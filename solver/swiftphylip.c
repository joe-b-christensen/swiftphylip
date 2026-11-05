/* SwiftPhylip: a byte-identical, single-pass reimplementation of vcf2phylip v2.9
   (PHYLIP, FASTA, NEXUS and binary NEXUS; deterministic modes).

   Output (files, stdout, exit status, and the temporary files vcf2phylip leaves behind when it
   raises) matches vcf2phylip.py v2.9 byte for byte on valid UTF-8 input.  The speed comes from
   reading the VCF once, turning genotypes into matrix characters with per-record lookup tables,
   and transposing in cache-sized tiles instead of re-reading a temporary file once per sample. */
#define _GNU_SOURCE
#include <ctype.h>
#include <dlfcn.h>
#include <errno.h>
#include <fcntl.h>
#include <limits.h>
#include <pthread.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/uio.h>
#include <unistd.h>
#include <zlib.h>
#if defined(__x86_64__)
#include <emmintrin.h>
#include <immintrin.h>
#endif

#define VERSION "1.0"

/* ------------------------------------------------------------------ utilities */
static void die(const char *msg) {
    fprintf(stderr, "swiftphylip: %s\n", msg);
    exit(1);
}
static void *xmalloc(size_t n) {
    void *p = malloc(n ? n : 1);
    if (!p) die("out of memory");
    return p;
}
static void *xrealloc(void *p, size_t n) {
    p = realloc(p, n ? n : 1);
    if (!p) die("out of memory");
    return p;
}

typedef struct { char *p; size_t n, cap; } Buf;
static inline void buf_grow(Buf *b, size_t need) {
    if (b->n + need > b->cap) {
        size_t c = b->cap ? b->cap : 4096;
        while (c < b->n + need) c *= 2;
        b->p = xrealloc(b->p, c);
        b->cap = c;
    }
}
static inline void buf_put(Buf *b, const void *s, size_t k) { buf_grow(b, k); memcpy(b->p + b->n, s, k); b->n += k; }
static inline void buf_puts(Buf *b, const char *s) { buf_put(b, s, strlen(s)); }
static inline void buf_putc(Buf *b, char c) { buf_grow(b, 1); b->p[b->n++] = c; }
static void buf_putll(Buf *b, long long v) {
    char t[32]; int k = 0; unsigned long long u = v < 0 ? (unsigned long long)(-(v + 1)) + 1 : (unsigned long long)v;
    do { t[k++] = (char)('0' + u % 10); u /= 10; } while (u);
    if (v < 0) buf_putc(b, '-');
    while (k) buf_putc(b, t[--k]);
}

static void write_all(int fd, const void *p, size_t n) {
    const char *c = p;
    while (n) {
        ssize_t w = write(fd, c, n);
        if (w < 0) { if (errno == EINTR) continue; die("write failed"); }
        c += w; n -= (size_t)w;
    }
}

/* stdout is buffered and written in large pieces; content equals Python's print output. */
static Buf OUT;
static void out_flush(void) { if (OUT.n) { write_all(1, OUT.p, OUT.n); OUT.n = 0; } }
static inline void out_maybe_flush(void) { if (OUT.n > (1u << 20)) out_flush(); }

/* Buffered file writer. */
typedef struct { int fd; Buf b; } FileW;
static void fw_open(FileW *w, const char *path) {
    memset(w, 0, sizeof *w);
    w->fd = open(path, O_WRONLY | O_CREAT | O_TRUNC, 0644);
    if (w->fd < 0) {
        out_flush();
        fprintf(stderr, "Traceback (most recent call last):\n  ...\nFileNotFoundError: [Errno %d] %s: '%s'\n", errno, strerror(errno), path);
        exit(1);
    }
}
static inline void fw_maybe_flush(FileW *w) { if (w->b.n > (1u << 20)) { write_all(w->fd, w->b.p, w->b.n); w->b.n = 0; } }
static void fw_close(FileW *w) { if (w->b.n) write_all(w->fd, w->b.p, w->b.n); w->b.n = 0; close(w->fd); w->fd = -1; }

/* ------------------------------------------------------------------ UTF-8 helpers */
static inline size_t cp_count(const unsigned char *p, size_t n) {
    size_t c = 0;
    for (size_t i = 0; i < n; i++) c += (p[i] & 0xC0) != 0x80;
    return c;
}
/* Python str.isspace() code points, encoded in UTF-8: length of the whitespace char at p (0 if none). */
static inline int ws_lead(const unsigned char *p, const unsigned char *e) {
    unsigned c = p[0];
    if (c == ' ' || (c >= 9 && c <= 13) || (c >= 0x1c && c <= 0x1f)) return 1;
    if (c < 0x80) return 0;
    if (c == 0xC2 && e - p >= 2 && (p[1] == 0x85 || p[1] == 0xA0)) return 2;
    if (e - p >= 3) {
        if (c == 0xE1 && p[1] == 0x9A && p[2] == 0x80) return 3;
        if (c == 0xE2 && p[1] == 0x80 && ((p[2] >= 0x80 && p[2] <= 0x8A) || p[2] == 0xA8 || p[2] == 0xA9 || p[2] == 0xAF)) return 3;
        if (c == 0xE2 && p[1] == 0x81 && p[2] == 0x9F) return 3;
        if (c == 0xE3 && p[1] == 0x80 && p[2] == 0x80) return 3;
    }
    return 0;
}
static inline int ws_trail(const unsigned char *b, const unsigned char *e) {
    unsigned c = e[-1];
    if (c == ' ' || (c >= 9 && c <= 13) || (c >= 0x1c && c <= 0x1f)) return 1;
    if (c < 0x80) return 0;
    if (e - b >= 2 && e[-2] == 0xC2 && (c == 0x85 || c == 0xA0)) return 2;
    if (e - b >= 3) {
        const unsigned char *p = e - 3;
        if (p[0] == 0xE1 && p[1] == 0x9A && p[2] == 0x80) return 3;
        if (p[0] == 0xE2 && p[1] == 0x80 && ((p[2] >= 0x80 && p[2] <= 0x8A) || p[2] == 0xA8 || p[2] == 0xA9 || p[2] == 0xAF)) return 3;
        if (p[0] == 0xE2 && p[1] == 0x81 && p[2] == 0x9F) return 3;
        if (p[0] == 0xE3 && p[1] == 0x80 && p[2] == 0x80) return 3;
    }
    return 0;
}

/* ------------------------------------------------------------------ input */
/* Decompressed input arrives as a sequence of chunks.  Plain files are memory-mapped (one chunk,
   no copy).  gzip input is decompressed with libdeflate when the shared library is present
   (loaded at run time; one member or a batch of BGZF members per chunk) and with zlib otherwise.
   With two or more threads, decompression runs in a producer thread ahead of the parser. */
typedef struct libdeflate_decompressor LDD;
static LDD *(*ld_alloc)(void);
static int (*ld_gzip_ex)(LDD *, const void *, size_t, void *, size_t, size_t *, size_t *);
static int ld_state = -1;
static int ld_available(void) {
    if (ld_state >= 0) return ld_state;
    ld_state = 0;
    const char *yes = getenv("SWIFTPHYLIP_LIBDEFLATE");   /* opt-in: slower than zlib on highly repetitive VCF text */
    if (!yes || !*yes || *yes == '0') return 0;
    void *h = dlopen("libdeflate.so.0", RTLD_NOW | RTLD_LOCAL);
    if (!h) return 0;
    ld_alloc = (LDD * (*)(void)) dlsym(h, "libdeflate_alloc_decompressor");
    ld_gzip_ex = (int (*)(LDD *, const void *, size_t, void *, size_t, size_t *, size_t *))dlsym(h, "libdeflate_gzip_decompress_ex");
    ld_state = ld_alloc && ld_gzip_ex;
    return ld_state;
}
static int THREADS = 1;
#include <time.h>
static int PROFILE;
static double T_START;
static double now_s(void) { struct timespec t; clock_gettime(CLOCK_MONOTONIC, &t); return t.tv_sec + t.tv_nsec * 1e-9; }
#include <sys/resource.h>
static void prof(const char *what) {
    if (!PROFILE) return;
    struct rusage ru; getrusage(RUSAGE_SELF, &ru);
    fprintf(stderr, "[prof] %-12s %8.2f ms  user %.1f sys %.1f minflt %ld\n", what, (now_s() - T_START) * 1e3,
            ru.ru_utime.tv_sec * 1e3 + ru.ru_utime.tv_usec / 1e3, ru.ru_stime.tv_sec * 1e3 + ru.ru_stime.tv_usec / 1e3, ru.ru_minflt);
}

static __thread int IN_PRODUCER;
static void gz_fail(const char *what) {
    if (!IN_PRODUCER) out_flush();
    fprintf(stderr, "Traceback (most recent call last):\n  ...\ngzip.BadGzipFile: %s\n", what);
    exit(1);
}
static void open_fail(const char *path) {
    out_flush();
    fprintf(stderr, "Traceback (most recent call last):\n  ...\nOSError: [Errno %d] %s: '%s'\n", errno, strerror(errno), path);
    exit(1);
}

typedef struct { char *p; size_t n; } Chunk;

/* whole-file view of the (compressed or plain) input */
typedef struct { const unsigned char *p; size_t n; int mapped; } FileView;
static void view_open(FileView *v, const char *path) {
    int fd = open(path, O_RDONLY);
    if (fd < 0) open_fail(path);
    struct stat st;
    v->mapped = 0;
    if (fstat(fd, &st) == 0 && S_ISREG(st.st_mode)) {
        v->n = (size_t)st.st_size;
        if (v->n == 0) { v->p = (const unsigned char *)""; close(fd); return; }
        void *m = mmap(NULL, v->n, PROT_READ, MAP_PRIVATE, fd, 0);
        if (m != MAP_FAILED) {
            madvise(m, v->n, MADV_SEQUENTIAL);
            v->p = m; v->mapped = 1; close(fd); return;
        }
    }
    Buf b = {0};
    for (;;) {
        buf_grow(&b, 1 << 20);
        ssize_t k = read(fd, b.p + b.n, b.cap - b.n);
        if (k < 0) { if (errno == EINTR) continue; die("read failed"); }
        if (k == 0) break;
        b.n += (size_t)k;
    }
    close(fd);
    v->p = (const unsigned char *)(b.p ? b.p : ""); v->n = b.n;
}
static void view_close(FileView *v) {
    if (v->mapped) munmap((void *)v->p, v->n);
    else if (v->n) free((void *)v->p);
}

/* Sequential gzip decoder over a FileView. */
typedef struct {
    FileView *v;
    size_t pos;          /* offset of the next unread compressed byte */
    int use_ld;
    LDD *ld;
    z_stream zs; int zactive;
    size_t chunk;        /* target chunk size */
    int done;
} GzDec;

static void gz_init(GzDec *g, FileView *v, size_t chunk, int allow_ld) {
    memset(g, 0, sizeof *g);
    g->v = v; g->chunk = chunk;
    g->use_ld = allow_ld && ld_available();
    if (g->use_ld) { g->ld = ld_alloc(); if (!g->ld) g->use_ld = 0; }
}
/* skip zero padding; returns 0 at end of input */
static int gz_member_start(GzDec *g) {
    while (g->pos < g->v->n && g->v->p[g->pos] == 0) g->pos++;
    if (g->pos >= g->v->n) return 0;
    if (g->v->p[g->pos] != 0x1f) gz_fail("Not a gzipped file");
    return 1;
}
/* BGZF block size of the member at pos, or 0 */
static size_t bgzf_bsize(const unsigned char *p, size_t n) {
    if (n < 18 || p[0] != 0x1f || p[1] != 0x8b || p[2] != 8 || !(p[3] & 4)) return 0;
    size_t xlen = p[10] | (p[11] << 8);
    if (12 + xlen > n) return 0;
    const unsigned char *x = p + 12, *xe = x + xlen;
    while (x + 4 <= xe) {
        size_t slen = x[2] | (x[3] << 8);
        if (x[0] == 'B' && x[1] == 'C' && slen == 2 && x + 6 <= xe) return (size_t)(x[4] | (x[5] << 8)) + 1;
        x += 4 + slen;
    }
    return 0;
}
/* zlib streaming into dst; returns bytes produced (0 at end) */
static size_t gz_zlib_read(GzDec *g, char *dst, size_t cap) {
    size_t got = 0;
    while (got < cap && !g->done) {
        if (!g->zactive) {
            if (!gz_member_start(g)) { g->done = 1; break; }
            memset(&g->zs, 0, sizeof g->zs);
            if (inflateInit2(&g->zs, 16 + 15) != Z_OK) die("zlib init failed");
            g->zactive = 1;
        }
        size_t avail = g->v->n - g->pos;
        if (avail == 0) gz_fail("Compressed file ended before the end-of-stream marker was reached");
        if (avail > (1u << 30)) avail = 1u << 30;
        g->zs.next_in = (unsigned char *)g->v->p + g->pos;
        g->zs.avail_in = (uInt)avail;
        size_t want = cap - got;
        if (want > (1u << 30)) want = 1u << 30;
        g->zs.next_out = (unsigned char *)dst + got;
        g->zs.avail_out = (uInt)want;
        int ret = inflate(&g->zs, Z_NO_FLUSH);
        g->pos += avail - g->zs.avail_in;
        got += want - g->zs.avail_out;
        if (ret == Z_STREAM_END) { inflateEnd(&g->zs); g->zactive = 0; }
        else if (ret == Z_BUF_ERROR && g->zs.avail_out > 0 && g->pos >= g->v->n) gz_fail("Compressed file ended before the end-of-stream marker was reached");
        else if (ret != Z_OK && ret != Z_BUF_ERROR) gz_fail("invalid compressed data");
    }
    return got;
}
/* Produce the next chunk into *buf (grown as needed). Returns bytes, 0 at end. */
static size_t gz_next(GzDec *g, Buf *buf) {
    buf->n = 0;
    if (!g->use_ld || g->zactive) {
        buf_grow(buf, g->chunk);
        size_t k = gz_zlib_read(g, buf->p, g->chunk);
        buf->n = k;
        return k;
    }
    while (!g->done) {
        if (!gz_member_start(g)) { g->done = 1; break; }
        const unsigned char *p = g->v->p + g->pos;
        size_t rem = g->v->n - g->pos;
        size_t bs = bgzf_bsize(p, rem);
        size_t outsz;
        if (bs && bs <= rem) outsz = (size_t)(p[bs - 4] | (p[bs - 3] << 8) | (p[bs - 2] << 16) | ((size_t)p[bs - 1] << 24));
        else {
            size_t isz = (size_t)(g->v->p[g->v->n - 4] | (g->v->p[g->v->n - 3] << 8) | (g->v->p[g->v->n - 2] << 16) | ((size_t)g->v->p[g->v->n - 1] << 24));
            outsz = rem < 18 ? 64 : isz;
            if (outsz < rem * 2) outsz = rem * 4 + 4096;
        }
        if (outsz > ((size_t)512 << 20)) { /* too large to hold: stream this member with zlib */
            buf_grow(buf, g->chunk);
            size_t k = gz_zlib_read(g, buf->p + buf->n, g->chunk);
            buf->n += k;
            return buf->n;
        }
        for (int attempt = 0;; attempt++) {
            buf_grow(buf, outsz + 1);
            size_t in_used = 0, out_got = 0;
            int r = ld_gzip_ex(g->ld, p, bs && bs <= rem ? bs : rem, buf->p + buf->n, outsz, &in_used, &out_got);
            if (r == 0) { g->pos += in_used; buf->n += out_got; break; }
            if (r == 3 && attempt < 6 && outsz < ((size_t)1 << 31)) { outsz *= 4; continue; }   /* insufficient space */
            if (r == 3) { buf_grow(buf, g->chunk); size_t k = gz_zlib_read(g, buf->p + buf->n, g->chunk); buf->n += k; return buf->n; }
            /* bad data: let zlib produce the same failure semantics (it decodes what it can first) */
            buf_grow(buf, g->chunk);
            size_t k = gz_zlib_read(g, buf->p + buf->n, g->chunk);
            buf->n += k;
            return buf->n;
        }
        if (buf->n >= g->chunk || (!bs && buf->n > 0)) break;   /* batch small BGZF members up to the chunk size */
    }
    return buf->n;
}
static void gz_close(GzDec *g) { if (g->zactive) inflateEnd(&g->zs); free(g->ld); }

/* Producer thread: a ring of decompressed buffers. */
#define NRING 4
typedef struct {
    GzDec g;
    Buf bufs[NRING];
    size_t lens[NRING];
    int head, tail, count, finished;
    pthread_mutex_t mu; pthread_cond_t cv;
    pthread_t th;
} Producer;
static void *producer_main(void *arg) {
    Producer *P = arg;
    IN_PRODUCER = 1;
    for (;;) {
        pthread_mutex_lock(&P->mu);
        while (P->count == NRING) pthread_cond_wait(&P->cv, &P->mu);
        int slot = P->tail;
        pthread_mutex_unlock(&P->mu);
        size_t k = gz_next(&P->g, &P->bufs[slot]);
        pthread_mutex_lock(&P->mu);
        P->lens[slot] = k;
        P->tail = (P->tail + 1) % NRING;
        P->count++;
        if (k == 0) P->finished = 1;
        pthread_cond_broadcast(&P->cv);
        pthread_mutex_unlock(&P->mu);
        if (k == 0) break;
    }
    return NULL;
}

typedef struct {
    FileView v;
    int gz, threaded;
    GzDec g; Buf gbuf;
    Producer *P; int held;   /* slot held by the consumer, -1 if none */
    int plain_done;
    /* line state */
    Chunk cur; size_t pos; int have, has_cr, pending_cr, eof;
    Buf carry; int carry_out;
    int borrowed;
} Lines;

static void lines_open(Lines *L, const char *path, int gz, size_t chunk, int header_pass) {
    memset(L, 0, sizeof *L);
    view_open(&L->v, path);
    L->gz = gz;
    L->held = -1;
    if (gz) {
        int allow_ld = !header_pass;
        if (!header_pass && THREADS >= 2) {
            L->P = calloc(1, sizeof(Producer));
            gz_init(&L->P->g, &L->v, chunk, allow_ld);
            pthread_mutex_init(&L->P->mu, NULL); pthread_cond_init(&L->P->cv, NULL);
            if (pthread_create(&L->P->th, NULL, producer_main, L->P) != 0) die("thread create failed");
            L->threaded = 1;
        } else gz_init(&L->g, &L->v, chunk, allow_ld);
    }
}
/* next raw chunk; 0 at end */
static int source_next(Lines *L, Chunk *c) {
    if (!L->gz) {
        if (L->plain_done) return 0;
        L->plain_done = 1;
        c->p = (char *)L->v.p; c->n = L->v.n;
        return c->n > 0;
    }
    if (!L->threaded) {
        size_t k = gz_next(&L->g, &L->gbuf);
        c->p = L->gbuf.p; c->n = k;
        return k > 0;
    }
    Producer *P = L->P;
    pthread_mutex_lock(&P->mu);
    if (L->held >= 0) { P->head = (P->head + 1) % NRING; P->count--; L->held = -1; pthread_cond_broadcast(&P->cv); }
    while (P->count == 0) pthread_cond_wait(&P->cv, &P->mu);
    int slot = P->head;
    size_t k = P->lens[slot];
    pthread_mutex_unlock(&P->mu);
    if (k == 0) return 0;
    L->held = slot;
    c->p = P->bufs[slot].p; c->n = k;
    return 1;
}

/* First '\n' or '\r' in [p, e), or NULL: one pass instead of memchr('\n') plus a '\r' scan. */
#if defined(__x86_64__)
__attribute__((target("avx2")))
static char *find_eol_avx2(char *p, char *e) {
    const __m256i nl = _mm256_set1_epi8('\n'), cr = _mm256_set1_epi8('\r');
    while (e - p >= 32) {
        __m256i v = _mm256_loadu_si256((const __m256i *)p);
        unsigned m = (unsigned)_mm256_movemask_epi8(_mm256_or_si256(_mm256_cmpeq_epi8(v, nl), _mm256_cmpeq_epi8(v, cr)));
        if (m) return p + __builtin_ctz(m);
        p += 32;
    }
    for (; p < e; p++) if (*p == '\n' || *p == '\r') return p;
    return NULL;
}
#endif
static int HAVE_AVX2_EOL;
static inline char *find_eol(char *p, char *e) {
#if defined(__x86_64__)
    if (HAVE_AVX2_EOL) return find_eol_avx2(p, e);
#endif
    char *nl = memchr(p, '\n', (size_t)(e - p));
    char *c = memchr(p, '\r', (size_t)((nl ? nl : e) - p));
    return c ? c : nl;
}
static inline int lines_next(Lines *L, char **ls, char **le) {
    if (L->carry_out) { L->carry.n = 0; L->carry_out = 0; }
    for (;;) {
        if (!L->have) {
            if (L->eof || !source_next(L, &L->cur)) {
                L->eof = 1;
                if (L->carry.n) { *ls = L->carry.p; *le = L->carry.p + L->carry.n; L->carry_out = 1; return 1; }
                return 0;
            }
            L->have = 1; L->pos = 0;
            if (L->pending_cr) { L->pending_cr = 0; if (L->cur.n && L->cur.p[0] == '\n') L->pos = 1; }
        }
        char *p = L->cur.p + L->pos, *e = L->cur.p + L->cur.n;
        char *t = find_eol(p, e);
        if (t) {
            size_t tl = 1;
            if (*t == '\r') { if (t + 1 < e) tl = t[1] == '\n' ? 2 : 1; else L->pending_cr = 1; }
            L->pos = (size_t)(t - L->cur.p) + tl;
            if (L->pos >= L->cur.n) L->have = 0;
            if (!L->carry.n) { *ls = p; *le = t; return 1; }
            buf_put(&L->carry, p, (size_t)(t - p));
            *ls = L->carry.p; *le = L->carry.p + L->carry.n; L->carry_out = 1;
            return 1;
        }
        buf_put(&L->carry, p, (size_t)(e - p));
        L->have = 0;
    }
}
static void lines_close(Lines *L) {
    if (L->threaded) {
        Producer *P = L->P;
        /* drain so the producer can finish */
        pthread_mutex_lock(&P->mu);
        while (!P->finished) {
            if (P->count) { P->head = (P->head + 1) % NRING; P->count--; pthread_cond_broadcast(&P->cv); }
            else pthread_cond_wait(&P->cv, &P->mu);
        }
        pthread_mutex_unlock(&P->mu);
        pthread_join(P->th, NULL);
        gz_close(&P->g);
        for (int i = 0; i < NRING; i++) free(P->bufs[i].p);
        free(P);
    } else if (L->gz) { gz_close(&L->g); free(L->gbuf.p); }
    free(L->carry.p);
    if (!L->borrowed) view_close(&L->v);
}
/* lines of an in-memory range (a part of a mapped plain file) */
static void lines_open_mem(Lines *L, const char *p, size_t n) {
    memset(L, 0, sizeof *L);
    L->v.p = (const unsigned char *)p; L->v.n = n;
    L->borrowed = 1;
    L->held = -1;
}

/* ------------------------------------------------------------------ IUPAC table (vcf2phylip AMBIG) */
static const char *AMBIG_KEYS[] = {
    "A","C","G","N","T","*A","*C","*G","*N","*T","AC","AG","AN","AT","CG","CN","CT","GN","GT","NT",
    "*AC","*AG","*AN","*AT","*CG","*CN","*CT","*GN","*GT","*NT","ACG","ACN","ACT","AGN","AGT","ANT","CGN","CGT","CNT","GNT",
    "*ACG","*ACN","*ACT","*AGN","*AGT","*ANT","*CGN","*CGT","*CNT","*GNT","ACGN","ACGT","ACNT","AGNT","CGNT",
    "*ACGN","*ACGT","*ACNT","*AGNT","*CGNT","*","*ACGNT"};
static const char AMBIG_VALS[] =
    "ACGNTacgntMRaWScYgKtmrawscygktVmHrDwsBykvmhrdwsbykvNhdbvNhdb-N";
static char AMBIG64[64];
static unsigned char AMBIGOK[64];
static int sym_bit(int c) {
    switch (c) { case '*': return 1; case 'A': return 2; case 'C': return 4; case 'G': return 8; case 'N': return 16; case 'T': return 32; }
    return 0;
}
static void ambig_init(void) {
    int nk = (int)(sizeof AMBIG_KEYS / sizeof AMBIG_KEYS[0]);
    for (int i = 0; i < nk; i++) {
        int m = 0;
        for (const char *k = AMBIG_KEYS[i]; *k; k++) m |= sym_bit(*k);
        AMBIG64[m] = AMBIG_VALS[i];
        AMBIGOK[m] = 1;
    }
}
static int ambig_lookup_str(const char *s, size_t n, char *out) {
    int nk = (int)(sizeof AMBIG_KEYS / sizeof AMBIG_KEYS[0]);
    for (int i = 0; i < nk; i++)
        if (strlen(AMBIG_KEYS[i]) == n && memcmp(AMBIG_KEYS[i], s, n) == 0) { *out = AMBIG_VALS[i]; return 1; }
    return 0;
}

/* ------------------------------------------------------------------ matrix storage */
/* Rows (one per accepted site, one byte per sample) are gathered RB at a time and transposed into
   sample-major tiles of B sites. Tiles stay in memory up to a budget, then go to an anonymous
   temporary file. */
#define RB 64
struct Matrix;
typedef struct { char *mem; off_t off; size_t nrec; struct Matrix *owner; } Tile;
typedef struct Matrix {
    int N;
    size_t B;
    char *rows;         /* RB x N row-major */
    int nrows;
    char *tile;         /* N x B sample-major, being filled */
    size_t in_tile;
    Tile *tiles; size_t ntiles, captiles;
    size_t total;
    size_t ram, budget;
    int spill_fd; off_t spill_off;
    const char *spill_dir;
    char *map;
} Matrix;


/* Tile buffers: 2 MB-aligned anonymous mappings advised for transparent huge pages, so a tile costs
   two huge-page faults instead of ~1000 small-page faults (first-touch faults dominate otherwise). */
static char *tile_alloc(size_t sz) {
    size_t al = 2u << 20;
    size_t len = (sz + 64 + al - 1) & ~(al - 1);
    char *p = mmap(NULL, len + al, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
    if (p == MAP_FAILED) return xmalloc(sz + 64);
    char *q = (char *)(((uintptr_t)p + al - 1) & ~(uintptr_t)(al - 1));
    if (q > p) munmap(p, (size_t)(q - p));
    size_t tail = (size_t)((p + len + al) - (q + len));
    if (tail) munmap(q + len, tail);
    madvise(q, len, MADV_HUGEPAGE);
    return q;
}
static void mat_init(Matrix *m, int N, size_t budget, const char *spill_dir) {
    memset(m, 0, sizeof *m);
    m->N = N;
    size_t B = ((2u << 20) - 64) / (size_t)(N > 0 ? N : 1);
    B = (B / RB) * RB;
    if (B < RB) B = RB;
    if (B > (1u << 16)) B = 1u << 16;
    m->B = B;
    m->rows = xmalloc((size_t)RB * N + 64);
    m->budget = budget;
    m->spill_fd = -1;
    m->spill_dir = spill_dir;
}
static inline char *mat_row(Matrix *m) { return m->rows + (size_t)m->nrows * m->N; }

static void transpose_rows(Matrix *m, int nrows) {
    const int N = m->N;
    const size_t B = m->B;
    char *dst0 = m->tile + m->in_tile;
    int s = 0;
#if defined(__x86_64__)
    if (nrows == RB) {
        for (; s + 16 <= N; s += 16) {
            for (int k0 = 0; k0 < RB; k0 += 16) {
                __m128i r[16];
                for (int i = 0; i < 16; i++) r[i] = _mm_loadu_si128((const __m128i *)(m->rows + (size_t)(k0 + i) * N + s));
                __m128i t[16];
                for (int i = 0; i < 8; i++) { t[2 * i] = _mm_unpacklo_epi8(r[2 * i], r[2 * i + 1]); t[2 * i + 1] = _mm_unpackhi_epi8(r[2 * i], r[2 * i + 1]); }
                __m128i u[16];
                for (int i = 0; i < 4; i++) {
                    u[4 * i + 0] = _mm_unpacklo_epi16(t[4 * i + 0], t[4 * i + 2]);
                    u[4 * i + 1] = _mm_unpackhi_epi16(t[4 * i + 0], t[4 * i + 2]);
                    u[4 * i + 2] = _mm_unpacklo_epi16(t[4 * i + 1], t[4 * i + 3]);
                    u[4 * i + 3] = _mm_unpackhi_epi16(t[4 * i + 1], t[4 * i + 3]);
                }
                __m128i v[16];
                for (int i = 0; i < 2; i++) {
                    for (int j = 0; j < 4; j++) {
                        v[8 * i + 2 * j] = _mm_unpacklo_epi32(u[8 * i + j], u[8 * i + 4 + j]);
                        v[8 * i + 2 * j + 1] = _mm_unpackhi_epi32(u[8 * i + j], u[8 * i + 4 + j]);
                    }
                }
                for (int j = 0; j < 8; j++) {
                    __m128i lo = _mm_unpacklo_epi64(v[j], v[8 + j]);
                    __m128i hi = _mm_unpackhi_epi64(v[j], v[8 + j]);
                    _mm_storeu_si128((__m128i *)(dst0 + (size_t)(s + 2 * j) * B + k0), lo);
                    _mm_storeu_si128((__m128i *)(dst0 + (size_t)(s + 2 * j + 1) * B + k0), hi);
                }
            }
        }
    }
#endif
    for (; s < N; s++) {
        char *d = dst0 + (size_t)s * B;
        for (int k = 0; k < nrows; k++) d[k] = m->rows[(size_t)k * N + s];
    }
}

static void mat_push_tile(Matrix *m) {
    if (m->ntiles == m->captiles) { m->captiles = m->captiles ? m->captiles * 2 : 64; m->tiles = xrealloc(m->tiles, m->captiles * sizeof(Tile)); }
    Tile *t = &m->tiles[m->ntiles++];
    t->nrec = m->in_tile;
    t->mem = m->tile;
    t->off = -1;
    t->owner = m;
    size_t sz = (size_t)m->N * m->B;
    m->ram += sz;
    if (m->ram > m->budget && m->in_tile == m->B) {
        if (m->spill_fd < 0) {
            m->spill_fd = open(m->spill_dir, O_TMPFILE | O_RDWR, 0600);
            if (m->spill_fd < 0) {
                char tmpl[4096];
                snprintf(tmpl, sizeof tmpl, "%s/.swiftphylip_spill_XXXXXX", m->spill_dir);
                m->spill_fd = mkstemp(tmpl);
                if (m->spill_fd < 0) die("cannot create spill file");
                unlink(tmpl);
            }
        }
        write_all(m->spill_fd, t->mem, sz);
        t->off = m->spill_off;
        m->spill_off += (off_t)sz;
        madvise(t->mem, sz, MADV_DONTNEED);   /* contents now live in the spill file */
        t->mem = NULL;
        m->ram -= sz;
    }
    m->tile = NULL;
    m->in_tile = 0;
}

static void mat_flush_rows(Matrix *m) {
    if (!m->nrows) return;
    if (!m->tile) m->tile = tile_alloc((size_t)m->N * m->B);
    transpose_rows(m, m->nrows);
    m->in_tile += (size_t)m->nrows;
    m->nrows = 0;
    if (m->in_tile == m->B) mat_push_tile(m);
}
static inline void mat_commit(Matrix *m) {
    m->total++;
    if (++m->nrows == RB) mat_flush_rows(m);
}
static void mat_finish(Matrix *m) {
    mat_flush_rows(m);
    if (m->in_tile) mat_push_tile(m);
}
/* Pointer to tile data (maps spilled tiles back). */
static const char *tile_data(Matrix *mm, Tile *t) {
    (void)mm;
    if (t->mem) return t->mem;
    Matrix *m = t->owner;
    if (!m->map) {
        m->map = mmap(NULL, (size_t)m->spill_off, PROT_READ, MAP_PRIVATE, m->spill_fd, 0);
        if (m->map == MAP_FAILED) die("mmap of spill file failed");
    }
    return m->map + t->off;
}
/* append the (finished) tiles of src after those of dst */
static void mat_append(Matrix *dst, Matrix *src) {
    for (size_t i = 0; i < src->ntiles; i++) {
        if (dst->ntiles == dst->captiles) { dst->captiles = dst->captiles ? dst->captiles * 2 : 64; dst->tiles = xrealloc(dst->tiles, dst->captiles * sizeof(Tile)); }
        dst->tiles[dst->ntiles++] = src->tiles[i];
    }
    dst->total += src->total;
}
/* Row-major dump of the matrix so far (pending rows included): vcf2phylip's .tmp layout. */
static void mat_dump_rowmajor(Matrix *m, const char *path) {
    FileW w; fw_open(&w, path);
    for (size_t ti = 0; ti < m->ntiles; ti++) {
        Tile *t = &m->tiles[ti];
        const char *d = tile_data(m, t);
        for (size_t k = 0; k < t->nrec; k++) {
            buf_grow(&w.b, (size_t)m->N + 1);
            for (int s = 0; s < m->N; s++) w.b.p[w.b.n++] = d[(size_t)s * m->B + k];
            w.b.p[w.b.n++] = '\n';
            fw_maybe_flush(&w);
        }
    }
    if (m->tile) {
        for (size_t k = 0; k < m->in_tile; k++) {
            buf_grow(&w.b, (size_t)m->N + 1);
            for (int s = 0; s < m->N; s++) w.b.p[w.b.n++] = m->tile[(size_t)s * m->B + k];
            w.b.p[w.b.n++] = '\n';
            fw_maybe_flush(&w);
        }
    }
    for (int k = 0; k < m->nrows; k++) {
        buf_put(&w.b, m->rows + (size_t)k * m->N, (size_t)m->N);
        buf_putc(&w.b, '\n');
        fw_maybe_flush(&w);
    }
    fw_close(&w);
}

/* ------------------------------------------------------------------ options (argparse emulation) */
typedef struct {
    const char *filename, *folder, *prefix, *outgroup;
    const char *msl_str;
    long long msl; int msl_set;
    int phylipdisable, fasta, nexus, nexusbin, write_used;
} Opts;

static void usage_error(const char *msg) {
    fprintf(stderr, "usage: vcf2phylip.py [-h] -i FILENAME [--output-folder FOLDER] [--output-prefix PREFIX] [-m MIN_SAMPLES_LOCUS] [-o OUTGROUP] [-p] [-f] [-n] [-b] [-r] [-w] [-v]\nvcf2phylip.py: error: %s\n", msg);
    exit(2);
}
/* Python int(): surrounding whitespace, optional sign, digits with single underscores between. */
static int parse_pyint(const char *s, long long *out) {
    const char *p = s;
    while (*p && strchr(" \t\n\r\f\v", *p)) p++;
    const char *e = p + strlen(p);
    while (e > p && strchr(" \t\n\r\f\v", e[-1])) e--;
    int neg = 0;
    if (p < e && (*p == '+' || *p == '-')) { neg = *p == '-'; p++; }
    if (p >= e || !isdigit((unsigned char)*p) || !isdigit((unsigned char)e[-1])) return 0;
    long long v = 0; int sat = 0;
    for (const char *q = p; q < e; q++) {
        if (*q == '_') { if (!isdigit((unsigned char)q[1])) return 0; continue; }
        if (!isdigit((unsigned char)*q)) return 0;
        if (v > (LLONG_MAX - 9) / 10) sat = 1; else v = v * 10 + (*q - '0');
    }
    if (sat) v = LLONG_MAX / 2;
    *out = neg ? -v : v;
    return 1;
}
static int looks_negative_number(const char *s) {
    /* argparse: '^-\d+$|^-\d*\.\d+$' */
    if (s[0] != '-') return 0;
    const char *p = s + 1;
    if (!*p) return 0;
    int digits = 0;
    while (isdigit((unsigned char)*p)) { p++; digits++; }
    if (!*p) return digits > 0;
    if (*p == '.') { p++; int d2 = 0; while (isdigit((unsigned char)*p)) { p++; d2++; } return !*p && d2 > 0; }
    return 0;
}

enum { O_INPUT, O_FOLDER, O_PREFIX, O_MSL, O_OUTGROUP, O_P, O_F, O_N, O_B, O_R, O_W, O_V, O_H, O_COUNT };
static const char *LONGS[O_COUNT] = {"--input", "--output-folder", "--output-prefix", "--min-samples-locus", "--outgroup",
                                     "--phylip-disable", "--fasta", "--nexus", "--nexus-binary", "--resolve-IUPAC",
                                     "--write-used-sites", "--version", "--help"};
static const char SHORTS[O_COUNT] = {'i', 0, 0, 'm', 'o', 'p', 'f', 'n', 'b', 'r', 'w', 'v', 'h'};
static int takes_value(int o) { return o == O_INPUT || o == O_FOLDER || o == O_PREFIX || o == O_MSL || o == O_OUTGROUP; }

static void apply_opt(Opts *op, int o, const char *val) {
    switch (o) {
    case O_INPUT: op->filename = val; break;
    case O_FOLDER: op->folder = val; break;
    case O_PREFIX: op->prefix = val; break;
    case O_OUTGROUP: op->outgroup = val; break;
    case O_MSL: {
        long long v;
        if (!parse_pyint(val, &v)) { char m[512]; snprintf(m, sizeof m, "argument -m/--min-samples-locus: invalid int value: '%s'", val); usage_error(m); }
        op->msl = v; op->msl_set = 1; op->msl_str = val; break;
    }
    case O_P: op->phylipdisable = 1; break;
    case O_F: op->fasta = 1; break;
    case O_N: op->nexus = 1; break;
    case O_B: op->nexusbin = 1; break;
    case O_W: op->write_used = 1; break;
    case O_R:
        fprintf(stderr, "swiftphylip: -r/--resolve-IUPAC (random resolution) is not supported; SwiftPhylip covers the deterministic modes only.\n");
        exit(2);
    case O_V: printf("vcf2phylip.py 2.9 (SwiftPhylip " VERSION ")\n"); exit(0);
    case O_H: printf("SwiftPhylip " VERSION ": byte-identical, fast vcf2phylip v2.9 (deterministic modes).\nusage: swiftphylip -i FILENAME [--output-folder FOLDER] [--output-prefix PREFIX] [-m MIN_SAMPLES_LOCUS] [-o OUTGROUP] [-p] [-f] [-n] [-b] [-w]\nEnvironment: SWIFTPHYLIP_THREADS (default: all CPUs), SWIFTPHYLIP_MEM_MB (tile memory before spilling, default 2048).\n"); exit(0);
    }
}

static void parse_args(int argc, char **argv, Opts *op) {
    memset(op, 0, sizeof *op);
    op->folder = "./";
    op->outgroup = "";
    op->msl = 4;
    for (int i = 1; i < argc; i++) {
        const char *a = argv[i];
        if (strcmp(a, "--") == 0) { if (i + 1 < argc) usage_error("unrecognized arguments"); break; }
        if (a[0] == '-' && a[1] == '-' && a[2]) {
            const char *eq = strchr(a, '=');
            size_t nlen = eq ? (size_t)(eq - a) : strlen(a);
            int found = -1, nmatch = 0;
            for (int o = 0; o < O_COUNT; o++) {
                if (strlen(LONGS[o]) == nlen && strncmp(LONGS[o], a, nlen) == 0) { found = o; nmatch = 1; break; }
            }
            if (found < 0) {
                for (int o = 0; o < O_COUNT; o++) if (strncmp(LONGS[o], a, nlen) == 0) { found = o; nmatch++; }
            }
            if (nmatch > 1) { char m[512]; snprintf(m, sizeof m, "ambiguous option: %.*s", (int)nlen, a); usage_error(m); }
            if (found < 0) { char m[512]; snprintf(m, sizeof m, "unrecognized arguments: %s", a); usage_error(m); }
            if (takes_value(found)) {
                const char *val;
                if (eq) val = eq + 1;
                else {
                    if (i + 1 >= argc || (argv[i + 1][0] == '-' && argv[i + 1][1] && !looks_negative_number(argv[i + 1]))) {
                        char m[512]; snprintf(m, sizeof m, "argument %s: expected one argument", LONGS[found]); usage_error(m);
                    }
                    val = argv[++i];
                }
                apply_opt(op, found, val);
            } else {
                if (eq) { char m[512]; snprintf(m, sizeof m, "argument %s: ignored explicit argument '%s'", LONGS[found], eq + 1); usage_error(m); }
                apply_opt(op, found, NULL);
            }
        } else if (a[0] == '-' && a[1] && a[1] != '-' && !looks_negative_number(a)) {
            for (const char *c = a + 1; *c; c++) {
                int found = -1;
                for (int o = 0; o < O_COUNT; o++) if (SHORTS[o] && SHORTS[o] == *c) found = o;
                if (found < 0) { char m[512]; snprintf(m, sizeof m, "unrecognized arguments: %s", a); usage_error(m); }
                if (takes_value(found)) {
                    const char *val;
                    if (c[1]) val = c + 1;
                    else {
                        if (i + 1 >= argc || (argv[i + 1][0] == '-' && argv[i + 1][1] && !looks_negative_number(argv[i + 1]))) {
                            char m[512]; snprintf(m, sizeof m, "argument -%c: expected one argument", *c); usage_error(m);
                        }
                        val = argv[++i];
                    }
                    apply_opt(op, found, val);
                    break;
                }
                apply_opt(op, found, NULL);
            }
        } else {
            char m[512]; snprintf(m, sizeof m, "unrecognized arguments: %s", a); usage_error(m);
        }
    }
    if (!op->filename) usage_error("the following arguments are required: -i/--input");
}

/* pathlib.PurePosixPath(folder, prefix) as a string */
static char *path_join_norm(const char *folder, const char *prefix) {
    Buf raw = {0};
    if (prefix[0] == '/') buf_puts(&raw, prefix);
    else { buf_puts(&raw, folder); buf_putc(&raw, '/'); buf_puts(&raw, prefix); }
    buf_putc(&raw, 0);
    const char *s = raw.p;
    Buf o = {0};
    size_t lead = 0;
    while (s[lead] == '/') lead++;
    if (lead == 2) buf_puts(&o, "//"); else if (lead >= 1) buf_putc(&o, '/');
    const char *p = s + lead;
    int first = 1;
    while (*p) {
        const char *q = strchr(p, '/');
        size_t len = q ? (size_t)(q - p) : strlen(p);
        if (len && !(len == 1 && p[0] == '.')) {
            if (!first) buf_putc(&o, '/');
            buf_put(&o, p, len);
            first = 0;
        }
        p += len;
        if (*p == '/') p++;
    }
    if (o.n == 0) buf_putc(&o, '.');
    buf_putc(&o, 0);
    free(raw.p);
    return o.p;
}
static void mkdir_parents(const char *path) {
    char *t = strdup(path);
    size_t n = strlen(t);
    for (size_t i = 1; i <= n; i++) {
        if (t[i] == '/' || t[i] == 0) {
            char c = t[i]; t[i] = 0;
            if (t[0]) mkdir(t, 0777);
            t[i] = c;
        }
    }
    free(t);
}

/* ------------------------------------------------------------------ main conversion */
typedef struct { const char *p; size_t n; } Str;

static Opts OPT;
static int N;
static Str *NAMES;
static size_t *NAMELEN_CP;
static char *OUTFILE;
static int NUC;          /* nucleotide matrices requested */
static Matrix NM, BM;
static FileW USED;
/* Per-thread parsing context. The main thread writes straight to stdout's buffer, the used-sites
   file and the global matrices; worker threads (parallel parsing of plain files) keep everything
   in memory and are merged in input order afterwards. */
typedef struct { long long li; size_t off, len; } Msg;
typedef struct {
    int is_main;
    Matrix *nm, *bm;
    Matrix own_nm, own_bm;
    Buf out; Msg *msgs; size_t nmsg, capmsg;
    Buf used;
    long long snp_num, snp_accepted, snp_shallow, mnp_num, snp_biallelic;
    int crashed; Buf crash_key;
    const char *p; size_t n;
} Ctx;
static Ctx MAINC;
static __thread Ctx *C;

static void outfile_path(Buf *b, const char *ext) { b->n = 0; buf_puts(b, OUTFILE); buf_puts(b, ext); buf_putc(b, 0); }

static void do_crash(const char *key, size_t keylen) {
    Buf pth = {0};
    if (NUC) { outfile_path(&pth, ".tmp"); mat_dump_rowmajor(&NM, pth.p); }
    if (OPT.nexusbin) { outfile_path(&pth, ".bin.tmp"); mat_dump_rowmajor(&BM, pth.p); }
    if (OPT.write_used) fw_close(&USED);
    out_flush();
    fprintf(stderr, "Traceback (most recent call last):\n  File \"vcf2phylip.py\", line 129, in get_matrix_column\n    column += AMBIG[geno_nuc]\nKeyError: '%.*s'\n", (int)keylen, key);
    exit(1);
}

static void print_skip(const char *b, const char *e) {
    if (C->is_main) {
        buf_puts(&OUT, "Skipping malformed line:\n");
        buf_put(&OUT, b, (size_t)(e - b));
        buf_putc(&OUT, '\n');
        out_maybe_flush();
        return;
    }
    if (C->nmsg == C->capmsg) { C->capmsg = C->capmsg ? 2 * C->capmsg : 64; C->msgs = xrealloc(C->msgs, C->capmsg * sizeof(Msg)); }
    Msg *m = &C->msgs[C->nmsg++];
    m->li = C->snp_num; m->off = C->out.n;
    buf_puts(&C->out, "Skipping malformed line:\n");
    buf_put(&C->out, b, (size_t)(e - b));
    buf_putc(&C->out, '\n');
    m->len = C->out.n - m->off;
}

/* per-record allele state */
#define C_BAD 0x80
#define C_POISON 0x40
static __thread unsigned char CLS[256];
static __thread Str ALLELE[11];          /* allele strings for digit ids 0..9 and '.'(=10) in the complex path */
static __thread char altbuf_static[256];

static void cls_reset(void) {
    memset(CLS, C_BAD, sizeof CLS);
    CLS['/'] = 0; CLS['|'] = 0; CLS['.'] = 16;
}

static int allele_class(const char *s, size_t n, int *complex) {
    if (n == 1) { int b = sym_bit((unsigned char)s[0]); if (b) return b; return C_POISON; }
    for (size_t i = 0; i < n; i++) if (!sym_bit((unsigned char)s[i])) return C_POISON;
    *complex = 1;                    /* empty or multi-character string of valid symbols */
    return 0;
}

static inline unsigned id_of(unsigned c) {
    if (c == '.') return 1u << 10;
    if (c >= '0' && c <= '9') return 1u << (c - '0');
    return 0;
}

static int str_cmp(const void *a, const void *b) {
    const Str *x = a, *y = b;
    size_t n = x->n < y->n ? x->n : y->n;
    int c = memcmp(x->p, y->p, n);
    if (c) return c;
    return (x->n > y->n) - (x->n < y->n);
}

/* Complex path: exact Python semantics on the sorted, de-duplicated allele strings.
   Returns 1 with *out set, 0 for KeyError (crash) with key in keybuf. */
static int complex_char(unsigned ids, char *out, Buf *keybuf) {
    Str v[11]; int nv = 0;
    for (int i = 0; i < 11; i++) if (ids & (1u << i)) v[nv++] = ALLELE[i];
    qsort(v, (size_t)nv, sizeof(Str), str_cmp);
    keybuf->n = 0;
    for (int i = 0; i < nv; i++) {
        if (i && v[i].n == v[i - 1].n && memcmp(v[i].p, v[i - 1].p, v[i].n) == 0) continue;
        buf_put(keybuf, v[i].p, v[i].n);
    }
    return ambig_lookup_str(keybuf->p ? keybuf->p : "", keybuf->n, out);
}


/* ------------------------------------------------------------------ AVX2 kernel for regular diploid lines */
/* Fields of the form "a?b\t" (one character, '/' or '|', one character, tab), 8 per 32-byte load.
   Each field becomes idx = 4*id(a) + id(b) with id('0')=0, id('1')=1, id('.')=2, other=3, and a
   16-entry per-record table gives the IUPAC character (0 = let the scalar code decide). Processing
   stops at the first group that is irregular or needs the scalar code; that code resumes there. */
static int HAVE_AVX2;
static __thread unsigned char T16[16];
static const unsigned char B16[16] = {'0', '1', '?', '?', '1', '2', '?', '?', '?', '?', '?', '?', '?', '?', '?', '?'};
#if defined(__x86_64__)
#include <immintrin.h>
__attribute__((target("avx2")))
static int simd_fields_avx2(const unsigned char *g, int N, int start, char *nrow, char *brow, long long *missing) {
    const int groups = (N - 1) / 8;
    const __m256i smask = _mm256_set1_epi32((int)0xFF00FF00u);
    const __m256i pat1 = _mm256_set1_epi32(0x09002F00), pat2 = _mm256_set1_epi32(0x09007C00);
    const __m256i c0 = _mm256_set1_epi8('0'), c1 = _mm256_set1_epi8('1'), cd = _mm256_set1_epi8('.');
    const __m256i one = _mm256_set1_epi8(1), two = _mm256_set1_epi8(2), three = _mm256_set1_epi8(3);
    const __m256i lo8 = _mm256_set1_epi32(0xFF), zero = _mm256_setzero_si256();
    const __m256i T = _mm256_broadcastsi128_si256(_mm_loadu_si128((const __m128i *)T16));
    const __m256i BT = _mm256_broadcastsi128_si256(_mm_loadu_si128((const __m128i *)B16));
    const __m256i pack = _mm256_setr_epi8(0, 4, 8, 12, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
                                          0, 4, 8, 12, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
    const __m256i perm = _mm256_setr_epi32(0, 4, 1, 1, 1, 1, 1, 1);
    long long miss = 0;
    int i;
    for (i = start / 8; i < groups; i++) {
        __m256i w = _mm256_loadu_si256((const __m256i *)(g + 32 * (size_t)i));
        __m256i m = _mm256_and_si256(w, smask);
        if (_mm256_movemask_epi8(_mm256_or_si256(_mm256_cmpeq_epi32(m, pat1), _mm256_cmpeq_epi32(m, pat2))) != -1) break;
        __m256i e1 = _mm256_cmpeq_epi8(w, c1), ed = _mm256_cmpeq_epi8(w, cd);
        __m256i any = _mm256_or_si256(_mm256_or_si256(_mm256_cmpeq_epi8(w, c0), e1), ed);
        __m256i id = _mm256_or_si256(_mm256_or_si256(_mm256_and_si256(e1, one), _mm256_and_si256(ed, two)), _mm256_andnot_si256(any, three));
        __m256i idx = _mm256_add_epi8(_mm256_slli_epi32(id, 2), _mm256_srli_epi32(id, 16));
        if (nrow) {
            __m256i ch = _mm256_shuffle_epi8(T, idx);
            if (_mm256_movemask_epi8(_mm256_cmpeq_epi32(_mm256_and_si256(ch, lo8), zero)) != 0) break;
            __m256i pk = _mm256_permutevar8x32_epi32(_mm256_shuffle_epi8(ch, pack), perm);
            _mm_storel_epi64((__m128i *)(nrow + 8 * (size_t)i), _mm256_castsi256_si128(pk));
        }
        if (brow) {
            __m256i pk = _mm256_permutevar8x32_epi32(_mm256_shuffle_epi8(_mm256_shuffle_epi8(BT, idx), pack), perm);
            _mm_storel_epi64((__m128i *)(brow + 8 * (size_t)i), _mm256_castsi256_si128(pk));
        }
        miss += __builtin_popcount((unsigned)_mm256_movemask_epi8(ed) & 0x11111111u);
    }
    *missing += miss;
    return i * 8 > start ? i * 8 : start;
}
static int HAVE_AVX512;
__attribute__((target("avx512f,avx512bw")))
static int simd_fields_avx512(const unsigned char *g, int N, char *nrow, char *brow, long long *missing) {
    const int groups = (N - 1) / 16;
    const __m512i smask = _mm512_set1_epi32((int)0xFF00FF00u);
    const __m512i pat1 = _mm512_set1_epi32(0x09002F00), pat2 = _mm512_set1_epi32(0x09007C00);
    const __m512i c0 = _mm512_set1_epi8('0'), c1 = _mm512_set1_epi8('1'), cd = _mm512_set1_epi8('.');
    const __m512i one = _mm512_set1_epi8(1), two = _mm512_set1_epi8(2), three = _mm512_set1_epi8(3), zero = _mm512_setzero_si512();
    const __m512i lo8 = _mm512_set1_epi32(0xFF);
    const __m512i T = _mm512_broadcast_i32x4(_mm_loadu_si128((const __m128i *)T16));
    const __m512i BT = _mm512_broadcast_i32x4(_mm_loadu_si128((const __m128i *)B16));
    long long miss = 0;
    int i;
    for (i = 0; i < groups; i++) {
        __m512i w = _mm512_loadu_si512((const void *)(g + 64 * (size_t)i));
        __m512i m = _mm512_and_si512(w, smask);
        if ((__mmask16)(_mm512_cmpeq_epi32_mask(m, pat1) | _mm512_cmpeq_epi32_mask(m, pat2)) != (__mmask16)0xFFFF) break;
        __mmask64 k0 = _mm512_cmpeq_epi8_mask(w, c0), k1 = _mm512_cmpeq_epi8_mask(w, c1), kd = _mm512_cmpeq_epi8_mask(w, cd);
        __m512i id = _mm512_mask_mov_epi8(three, k0, zero);
        id = _mm512_mask_mov_epi8(id, k1, one);
        id = _mm512_mask_mov_epi8(id, kd, two);
        __m512i idx = _mm512_add_epi8(_mm512_slli_epi32(id, 2), _mm512_srli_epi32(id, 16));
        if (nrow) {
            __m512i ch = _mm512_shuffle_epi8(T, idx);
            if (_mm512_cmpeq_epi32_mask(_mm512_and_si512(ch, lo8), zero)) break;
            _mm_storeu_si128((__m128i *)(nrow + 16 * (size_t)i), _mm512_cvtepi32_epi8(ch));
        }
        if (brow) _mm_storeu_si128((__m128i *)(brow + 16 * (size_t)i), _mm512_cvtepi32_epi8(_mm512_shuffle_epi8(BT, idx)));
        miss += __builtin_popcountll(kd & 0x1111111111111111ull);
    }
    *missing += miss;
    return i * 16;
}
#endif
static void build_T16(void) {
    static const unsigned char chs[3] = {'0', '1', '.'};
    for (int a = 0; a < 4; a++)
        for (int b = 0; b < 4; b++) {
            unsigned char r = 0;
            if (a < 3 && b < 3) {
                unsigned v = CLS[chs[a]] | CLS[chs[b]];
                if (!(v & (C_BAD | C_POISON)) && AMBIGOK[v & 63]) r = (unsigned char)AMBIG64[v & 63];
            }
            T16[4 * a + b] = r;
        }
}


/* f[1..9] = starts of fields 1..9 (just after the first nine tabs); 0 if there are fewer tabs. */
#if defined(__x86_64__)
__attribute__((target("avx2")))
static int find_tabs9_avx2(const unsigned char *p, const unsigned char *e, const unsigned char **f) {
    const __m256i tab = _mm256_set1_epi8('\t');
    int k = 1;
    const unsigned char *q = p;
    while (e - q >= 32) {
        unsigned m = (unsigned)_mm256_movemask_epi8(_mm256_cmpeq_epi8(_mm256_loadu_si256((const __m256i *)q), tab));
        while (m) { f[k++] = q + __builtin_ctz(m) + 1; if (k == 10) return 1; m &= m - 1; }
        q += 32;
    }
    for (; q < e; q++) if (*q == '\t') { f[k++] = q + 1; if (k == 10) return 1; }
    return 0;
}
#endif
static int find_tabs9(const unsigned char *p, const unsigned char *e, const unsigned char **f) {
#if defined(__x86_64__)
    if (HAVE_AVX2_EOL) return find_tabs9_avx2(p, e, f);
#endif
    for (int i = 1; i <= 9; i++) {
        const unsigned char *t = memchr(p, '\t', (size_t)(e - p));
        if (!t) return 0;
        f[i] = t + 1;
        p = t + 1;
    }
    return 1;
}
static int process_line(char *ls, char *le) {
    const unsigned char *b = (const unsigned char *)ls, *e = (const unsigned char *)le;
    int k;
    while (b < e && (k = ws_lead(b, e))) b += k;
    while (e > b && (k = ws_trail(b, e))) e -= k;
    if (b == e || *b == '#') return 0;
    C->snp_num++;
    if (C->is_main && C->snp_num % 500000 == 0) { buf_putll(&OUT, C->snp_num); buf_puts(&OUT, " genotypes processed.\n"); }
    /* first nine fields */
    const unsigned char *f[10];
    f[0] = b;
    const unsigned char *p = b;
    if (!find_tabs9(b, e, f)) { print_skip((const char *)b, (const char *)e); return 0; }   /* fewer than 10 fields */
    const unsigned char *ref = f[3], *refe = f[4] - 1, *alt = f[4], *alte = f[5] - 1;
    size_t refn = (size_t)(refe - ref), altn = (size_t)(alte - alt);
    /* is_snp */
    int snp = 0;
    size_t alt_cp = 0;
    int simple = refn == 1 && altn == 1 && ref[0] < 0x80 && alt[0] < 0x80 && alt[0] != ',';
    size_t ref_cp = simple ? 1 : cp_count(ref, refn);
    if (simple) { snp = 1; alt_cp = 1; }
    else if (ref_cp == 1) {
        size_t commas = 0, nonref = 0;
        alt_cp = cp_count(alt, altn);
        for (size_t i = 0; i < altn; i++) commas += alt[i] == ',';
        if (altn >= 9) {
            const unsigned char *q = alt;
            while ((q = memmem(q, (size_t)(alte - q), "<NON_REF>", 9))) { nonref++; q += 9; }
        }
        size_t L = alt_cp - 9 * nonref + nonref;
        size_t c = commas + (refn == 1 && ref[0] == ',' ? nonref : 0);
        snp = (L - c == c + 1);
    }
    const unsigned char *g = f[9];
    int bin_rec = snp && OPT.nexusbin && alt_cp == 1;
    int nuc_rec = snp && NUC;
    int complex = 0;
    if (nuc_rec && simple) {
        /* one ASCII REF and one ASCII ALT character: classes for '0' and '1' only */
        unsigned char r0 = ref[0], a1 = alt[0];
        r0 = r0 == '-' ? '*' : (unsigned char)toupper(r0);
        a1 = a1 == '-' ? '*' : a1;
        int b0 = sym_bit(r0), b1 = sym_bit(a1);
        memset(CLS + '2', C_BAD, 8);
        CLS['0'] = (unsigned char)(b0 ? b0 : C_POISON);
        CLS['1'] = (unsigned char)(b1 ? b1 : C_POISON);
    } else if (nuc_rec) {
        memset(CLS + '0', C_BAD, 10);
        CLS['.'] = 16;
        /* allele 0: REF.replace("-","*").upper(); non-ASCII -> poison */
        char nt0;
        int nt0_poison = 0;
        if (ref[0] < 0x80) { nt0 = (char)ref[0]; if (nt0 == '-') nt0 = '*'; else nt0 = (char)toupper((unsigned char)nt0); }
        else { nt0 = (char)0xFF; nt0_poison = 1; }
        (void)nt0_poison;
        /* ALT.replace("-","*").replace("<NON_REF>", nt0).split(",") */
        char *ab = altbuf_static;
        Buf dyn = {0};
        if (altn + 1 > sizeof altbuf_static) { buf_grow(&dyn, altn + 1); ab = dyn.p; }
        size_t an = 0;
        for (size_t i = 0; i < altn;) {
            if (alt[i] == '<' && altn - i >= 9 && memcmp(alt + i, "<NON_REF>", 9) == 0) { ab[an++] = nt0; i += 9; }
            else { ab[an++] = alt[i] == '-' ? '*' : (char)alt[i]; i++; }
        }
        static __thread char nt0buf[1];
        nt0buf[0] = nt0;
        ALLELE[0].p = nt0buf; ALLELE[0].n = 1;
        ALLELE[10].p = "N"; ALLELE[10].n = 1;
        CLS['0'] = (unsigned char)allele_class(nt0buf, 1, &complex);
        size_t st = 0; int idx = 1;
        for (size_t i = 0; i <= an; i++) {
            if (i == an || ab[i] == ',') {
                if (idx <= 9) {
                    ALLELE[idx].p = ab + st; ALLELE[idx].n = i - st;
                    CLS['0' + idx] = (unsigned char)allele_class(ab + st, i - st, &complex);
                }
                idx++;
                st = i + 1;
            }
        }
        if (complex) {
            /* copy allele strings out of the ALT buffer before it is reused */
            static __thread Buf keep;
            keep.n = 0;
            for (int i = 1; i <= 9 && i < idx; i++) buf_put(&keep, ALLELE[i].p, ALLELE[i].n);
            size_t o = 0;
            for (int i = 1; i <= 9 && i < idx; i++) { ALLELE[i].p = keep.p + o; o += ALLELE[i].n; }
        }
        free(dyn.p);
    }

    int s = 0;
    long long missing = 0;
    if (!nuc_rec && !bin_rec) {
        /* only field count and missing count are needed */
        p = g;
        for (;;) {
            if (p < e && *p == '.') missing++;
            const unsigned char *t = memchr(p, '\t', (size_t)(e - p));
            s++;
            if (!t) break;
            p = t + 1;
        }
        if (s != N) { print_skip((const char *)b, (const char *)e); return 0; }
        if (N - missing < OPT.msl) { C->snp_shallow++; return 0; }
        if (!snp) { C->mnp_num++; return 0; }
        return 0;   /* SNP, but neither nucleotide nor binary matrix requested */
    }

    char *nrow = nuc_rec ? mat_row(C->nm) : NULL;
    char *brow = bin_rec ? mat_row(C->bm) : NULL;
    int event = 0;                      /* 1 malformed, 2 crash */
    int ev_sample = -1;
    unsigned ev_ids = 0; unsigned ev_v = 0;
    p = g;
#if defined(__x86_64__)
    if (HAVE_AVX2 && !complex && N >= 9 && (size_t)(e - g) == 4 * (size_t)N - 1) {
        if (nrow) build_T16();
        int done = HAVE_AVX512 ? simd_fields_avx512(g, N, nrow, brow, &missing) : 0;
        done = simd_fields_avx2(g, N, done, nrow, brow, &missing);
        s = done;
        p = g + 4 * (size_t)done;
    }
#endif
    for (;;) {
        const unsigned char *q;
        if (p < e && *p == '.') missing++;
        if (p + 3 <= e && (p + 3 == e || p[3] == '\t') && (p[1] == '/' || p[1] == '|')) {
            q = p + 3;
            if (s < N) {
                unsigned c0 = p[0], c2 = p[2];
                if (nrow && !event) {
                    unsigned v = CLS[c0] | CLS[c2];
                    if (complex) {
                        if (v & C_BAD) { event = 1; ev_sample = s; }
                        else {
                            unsigned ids = id_of(c0) | id_of(c2);
                            char ch; static __thread Buf kb;
                            if (complex_char(ids, &ch, &kb)) nrow[s] = ch; else { event = 2; ev_sample = s; ev_ids = ids; }
                        }
                    } else if (v & C_BAD) { event = 1; ev_sample = s; }
                    else if ((v & C_POISON) || !AMBIGOK[v & 63]) { event = 2; ev_sample = s; ev_v = v; ev_ids = id_of(c0) | id_of(c2); }
                    else nrow[s] = AMBIG64[v & 63];
                }
                if (brow) brow[s] = ((c0 == '0' || c0 == '1') && (c2 == '0' || c2 == '1')) ? (char)('0' + (c0 - '0') + (c2 - '0')) : '?';
            }
        } else {
            const unsigned char *r = p;
            unsigned v = 0, ids = 0;
            while (r < e && *r != '\t' && *r != ':') {
                v |= CLS[*r];
                ids |= id_of(*r);
                r++;
            }
            size_t glen = (size_t)(r - p);
            if (r < e && *r == ':') { const unsigned char *t = memchr(r, '\t', (size_t)(e - r)); r = t ? t : e; }
            q = r;
            if (s < N) {
                if (nrow && !event) {
                    if (v & C_BAD) { event = 1; ev_sample = s; }
                    else if (complex) {
                        char ch; static __thread Buf kb;
                        if (complex_char(ids, &ch, &kb)) nrow[s] = ch; else { event = 2; ev_sample = s; ev_ids = ids; }
                    } else if ((v & C_POISON) || !AMBIGOK[v & 63]) { event = 2; ev_sample = s; ev_v = v; ev_ids = ids; }
                    else nrow[s] = AMBIG64[v & 63];
                }
                if (brow) {
                    char bc = '?';
                    if (glen == 3 && (p[1] == '/' || p[1] == '|') && (p[0] == '0' || p[0] == '1') && (p[2] == '0' || p[2] == '1'))
                        bc = (char)('0' + (p[0] - '0') + (p[2] - '0'));
                    brow[s] = bc;
                }
            }
        }
        s++;
        if (q >= e) break;
        p = q + 1;
    }
    (void)ev_sample; (void)ev_v;
    if (s != N) { print_skip((const char *)b, (const char *)e); return 0; }
    long long nsl = N - missing;
    if (nsl < OPT.msl) { C->snp_shallow++; return 0; }
    if (nuc_rec) {
        if (event == 1) { print_skip((const char *)b, (const char *)e); return 0; }
        if (event == 2) {
            /* the key Python reports (stderr only; not part of the compared output) */
            static __thread Buf kb;
            char ch;
            if (complex) complex_char(ev_ids, &ch, &kb);
            else {
                kb.n = 0;
                static const char sy[] = "*ACGNT";
                for (int i = 0; i < 6; i++) if (ev_v & (1u << i)) buf_putc(&kb, sy[i]);
            }
            if (C->is_main) do_crash(kb.p ? kb.p : "", kb.n);
            C->crashed = 1;
            C->crash_key.n = 0;
            buf_put(&C->crash_key, kb.p ? kb.p : "", kb.n);
            return 1;
        }
        C->snp_accepted++;
        mat_commit(C->nm);
        if (OPT.write_used) {
            Buf *ub = C->is_main ? &USED.b : &C->used;
            buf_put(ub, f[0], (size_t)(f[1] - 1 - f[0]));
            buf_putc(ub, '\t');
            buf_put(ub, f[1], (size_t)(f[2] - 1 - f[1]));
            buf_putc(ub, '\t');
            buf_putll(ub, nsl);
            buf_putc(ub, '\n');
            if (C->is_main) fw_maybe_flush(&USED);
        }
    }
    if (bin_rec) { C->snp_biallelic++; mat_commit(C->bm); }
    return 0;
}

/* ------------------------------------------------------------------ output */
static void write_pad(Buf *b, int s, size_t longest) {
    size_t pad = longest + 3 - NAMELEN_CP[s];
    buf_put(b, NAMES[s].p, NAMES[s].n);
    buf_grow(b, pad);
    memset(b->p + b->n, ' ', pad);
    b->n += pad;
}

/* Zero-copy output: each row is written with writev straight from the tiles. */
typedef struct { int fd; struct iovec iov[1024]; int n; size_t bytes; } IovW;
static void iov_flush(IovW *w) {
    struct iovec *v = w->iov;
    int n = w->n;
    while (n > 0) {
        ssize_t k = writev(w->fd, v, n > 1024 ? 1024 : n);
        if (k < 0) { if (errno == EINTR) continue; die("write failed"); }
        size_t kk = (size_t)k;
        while (n > 0 && kk >= v->iov_len) { kk -= v->iov_len; v++; n--; }
        if (n > 0 && kk) { v->iov_base = (char *)v->iov_base + kk; v->iov_len -= kk; }
    }
    w->n = 0;
}
static inline void iov_add(IovW *w, const void *p, size_t len) {
    if (!len) return;
    if (w->n == 1024) iov_flush(w);
    w->iov[w->n].iov_base = (void *)p;
    w->iov[w->n].iov_len = len;
    w->n++;
}
typedef struct {
    const char *path; Matrix *m; int fasta;
    Buf head; const char *tail;
} OutJob;
static int *ORDER; static int NORDER;
static Buf PREFIX_PAD, PREFIX_FA; static size_t *POFF_PAD, *POFF_FA;
static void run_job(OutJob *j) {
    IovW *w = xmalloc(sizeof(IovW));
    w->fd = open(j->path, O_WRONLY | O_CREAT | O_TRUNC, 0644);
    w->n = 0;
    if (w->fd < 0) die("cannot create output file");
    static const char nl = '\n';
    iov_add(w, j->head.p, j->head.n);
    Matrix *m = j->m;
    const char **td = xmalloc(sizeof(char *) * (m->ntiles + 1));
    for (size_t t = 0; t < m->ntiles; t++) td[t] = tile_data(m, &m->tiles[t]);
    for (int i = 0; i < NORDER; i++) {
        int s = ORDER[i];
        if (j->fasta) iov_add(w, PREFIX_FA.p + POFF_FA[s], POFF_FA[s + 1] - POFF_FA[s]);
        else iov_add(w, PREFIX_PAD.p + POFF_PAD[s], POFF_PAD[s + 1] - POFF_PAD[s]);
        for (size_t t = 0; t < m->ntiles; t++) iov_add(w, td[t] + (size_t)s * m->B, m->tiles[t].nrec);
        iov_add(w, &nl, 1);
    }
    if (j->tail) iov_add(w, j->tail, strlen(j->tail));
    iov_flush(w);
    close(w->fd);
    free(td); free(w);
}

static void write_matrices(int idx_out, size_t longest) {
    ORDER = xmalloc(sizeof(int) * (size_t)N);
    NORDER = 0;
    if (idx_out >= 0) ORDER[NORDER++] = idx_out;
    for (int s = 0; s < N; s++) if (s != idx_out) ORDER[NORDER++] = s;
    POFF_PAD = xmalloc(sizeof(size_t) * (size_t)(N + 1));
    POFF_FA = xmalloc(sizeof(size_t) * (size_t)(N + 1));
    for (int s = 0; s < N; s++) {
        POFF_PAD[s] = PREFIX_PAD.n; write_pad(&PREFIX_PAD, s, longest);
        POFF_FA[s] = PREFIX_FA.n; buf_putc(&PREFIX_FA, '>'); buf_put(&PREFIX_FA, NAMES[s].p, NAMES[s].n); buf_putc(&PREFIX_FA, '\n');
    }
    POFF_PAD[N] = PREFIX_PAD.n; POFF_FA[N] = PREFIX_FA.n;
    OutJob jobs[4]; int nj = 0;
    Buf pth = {0};
    memset(jobs, 0, sizeof jobs);
    if (!OPT.phylipdisable) {
        OutJob *j = &jobs[nj++]; outfile_path(&pth, ".phy"); j->path = strdup(pth.p); j->m = &NM;
        buf_putll(&j->head, N); buf_putc(&j->head, ' '); buf_putll(&j->head, MAINC.snp_accepted); buf_putc(&j->head, '\n');
    }
    if (OPT.fasta) { OutJob *j = &jobs[nj++]; outfile_path(&pth, ".fasta"); j->path = strdup(pth.p); j->m = &NM; j->fasta = 1; }
    if (OPT.nexus) {
        OutJob *j = &jobs[nj++]; outfile_path(&pth, ".nexus"); j->path = strdup(pth.p); j->m = &NM; j->tail = ";\nEND;\n";
        buf_puts(&j->head, "#NEXUS\n\nBEGIN DATA;\n\tDIMENSIONS NTAX="); buf_putll(&j->head, N);
        buf_puts(&j->head, " NCHAR="); buf_putll(&j->head, MAINC.snp_accepted);
        buf_puts(&j->head, ";\n\tFORMAT DATATYPE=DNA MISSING=N GAP=- ;\nMATRIX\n");
    }
    if (OPT.nexusbin) {
        OutJob *j = &jobs[nj++]; outfile_path(&pth, ".bin.nexus"); j->path = strdup(pth.p); j->m = &BM; j->tail = ";\nEND;\n";
        buf_puts(&j->head, "#NEXUS\n\nBEGIN DATA;\n\tDIMENSIONS NTAX="); buf_putll(&j->head, N);
        buf_puts(&j->head, " NCHAR="); buf_putll(&j->head, MAINC.snp_biallelic);
        buf_puts(&j->head, ";\n\tFORMAT DATATYPE=SNP MISSING=? GAP=- ;\nMATRIX\n");
    }
    /* map spilled tiles before any job runs */
    for (int i = 0; i < nj; i++) for (size_t t = 0; t < jobs[i].m->ntiles; t++) tile_data(jobs[i].m, &jobs[i].m->tiles[t]);
    for (int i = 0; i < nj; i++) run_job(&jobs[i]);
    free(pth.p);
}


/* ------------------------------------------------------------------ parallel parsing of plain files */
/* The mapped file is cut after a newline near each 1/T point (always a line boundary under
   universal newlines). Part 0 is parsed by the main thread exactly as in the serial path; parts
   1..T-1 are parsed into private contexts and merged in order: their messages are replayed with
   the progress lines re-inserted at the global counts, their used-site rows appended, their tiles
   appended, and a KeyError in a part is raised after everything before it has been merged. */
static void *part_main(void *arg) {
    Ctx *c = arg;
    C = c;
    cls_reset();
    Lines L; lines_open_mem(&L, c->p, c->n);
    char *ls, *le;
    while (lines_next(&L, &ls, &le)) if (process_line(ls, le)) break;
    lines_close(&L);
    return NULL;
}
static void parse_plain_parallel(const char *fn, size_t budget) {
    FileView v; view_open(&v, fn);
    int T = THREADS > 16 ? 16 : THREADS;
    if (v.n < (4u << 20)) T = 1;
    size_t cut[17];
    cut[0] = 0;
    for (int k = 1; k < T; k++) {
        size_t at = v.n * (size_t)k / (size_t)T;
        if (at < cut[k - 1]) at = cut[k - 1];
        const char *nl = at < v.n ? memchr(v.p + at, '\n', v.n - at) : NULL;
        cut[k] = nl ? (size_t)(nl - (const char *)v.p) + 1 : v.n;
    }
    cut[T] = v.n;
    Ctx *ctx = calloc((size_t)T, sizeof(Ctx));
    pthread_t th[17];
    for (int k = 1; k < T; k++) {
        Ctx *c = &ctx[k];
        c->p = (const char *)v.p + cut[k]; c->n = cut[k + 1] - cut[k];
        c->nm = &c->own_nm; c->bm = &c->own_bm;
        if (NUC) mat_init(c->nm, N, (OPT.nexusbin ? budget / 2 : budget) / (size_t)T, NM.spill_dir);
        if (OPT.nexusbin) mat_init(c->bm, N, (NUC ? budget / 2 : budget) / (size_t)T, BM.spill_dir);
        if (pthread_create(&th[k], NULL, part_main, c) != 0) die("thread create failed");
    }
    {
        Lines L; lines_open_mem(&L, (const char *)v.p, cut[1]);
        char *ls, *le;
        while (lines_next(&L, &ls, &le)) process_line(ls, le);
        lines_close(&L);
    }
    prof("parse0");
    if (NUC) mat_finish(&NM);
    if (OPT.nexusbin) mat_finish(&BM);
    for (int k = 1; k < T; k++) {
        pthread_join(th[k], NULL);
        Ctx *c = &ctx[k];
        long long base = MAINC.snp_num;
        long long next = (base / 500000 + 1) * 500000;
        for (size_t i = 0; i < c->nmsg; i++) {
            while (next <= base + c->msgs[i].li) { buf_putll(&OUT, next); buf_puts(&OUT, " genotypes processed.\n"); next += 500000; }
            buf_put(&OUT, c->out.p + c->msgs[i].off, c->msgs[i].len);
            out_maybe_flush();
        }
        while (next <= base + c->snp_num) { buf_putll(&OUT, next); buf_puts(&OUT, " genotypes processed.\n"); next += 500000; }
        if (OPT.write_used && c->used.n) { buf_put(&USED.b, c->used.p, c->used.n); fw_maybe_flush(&USED); }
        MAINC.snp_num += c->snp_num; MAINC.snp_accepted += c->snp_accepted; MAINC.snp_shallow += c->snp_shallow;
        MAINC.mnp_num += c->mnp_num; MAINC.snp_biallelic += c->snp_biallelic;
        if (NUC) { mat_finish(c->nm); mat_append(&NM, c->nm); }
        if (OPT.nexusbin) { mat_finish(c->bm); mat_append(&BM, c->bm); }
        if (c->crashed) do_crash(c->crash_key.p ? c->crash_key.p : "", c->crash_key.n);
    }
    prof("parse");
    /* the mapping stays alive: rows were already transposed into tiles, so nothing points into it */
    view_close(&v);
}

int main(int argc, char **argv) {
    T_START = now_s();
    PROFILE = getenv("SWIFTPHYLIP_PROFILE") != NULL;
    parse_args(argc, argv, &OPT);
    {
        long nc = sysconf(_SC_NPROCESSORS_ONLN);
        THREADS = nc > 0 ? (int)nc : 1;
        const char *t = getenv("SWIFTPHYLIP_THREADS");
        if (t && atoi(t) > 0) THREADS = atoi(t);
    }
    ambig_init();
    cls_reset();
#if defined(__x86_64__)
    HAVE_AVX2 = __builtin_cpu_supports("avx2");
    HAVE_AVX2_EOL = HAVE_AVX2;
    HAVE_AVX512 = HAVE_AVX2 && __builtin_cpu_supports("avx512f") && __builtin_cpu_supports("avx512bw");
    { const char *na = getenv("SWIFTPHYLIP_NO_AVX512"); if (na && *na && *na != '0') HAVE_AVX512 = 0; }
    { const char *na = getenv("SWIFTPHYLIP_NO_SIMD"); if (na && *na && *na != '0') HAVE_AVX2 = 0; }
#endif
    const char *fn = OPT.filename;
    size_t fnl = strlen(fn);
    int gz = fnl >= 3 && (fn[fnl - 3] == '.') && tolower((unsigned char)fn[fnl - 2]) == 'g' && tolower((unsigned char)fn[fnl - 1]) == 'z';
    char *og = strdup(OPT.outgroup);
    og[strcspn(og, ",;")] = 0;

    struct stat stt;
    if (stat(fn, &stt) != 0) {
        buf_puts(&OUT, "\nInput VCF file not found, please verify the provided path\n");
        out_flush();
        return 0;
    }
    /* pass 1: sample names from the first line starting with #CHROM */
    Buf hdr = {0};
    {
        Lines L; lines_open(&L, fn, gz, 64u << 10, 1);
        char *ls, *le;
        int found = 0;
        while (lines_next(&L, &ls, &le)) {
            if (le - ls >= 6 && memcmp(ls, "#CHROM", 6) == 0) { buf_put(&hdr, ls, (size_t)(le - ls)); found = 1; break; }
        }
        lines_close(&L);
        (void)found;
    }
    /* split the header on tabs; names = fields 9.. with every "./" removed */
    {
        size_t nf = 0;
        size_t cap = 16;
        Str *fields = xmalloc(cap * sizeof(Str));
        const char *p = hdr.p ? hdr.p : "", *e = p + hdr.n;
        for (;;) {
            const char *t = memchr(p, '\t', (size_t)(e - p));
            if (nf == cap) { cap *= 2; fields = xrealloc(fields, cap * sizeof(Str)); }
            fields[nf].p = p; fields[nf].n = t ? (size_t)(t - p) : (size_t)(e - p); nf++;
            if (!t) break;
            p = t + 1;
        }
        N = nf > 9 ? (int)(nf - 9) : 0;
        if (!hdr.n) N = 0;
        NAMES = xmalloc(sizeof(Str) * (size_t)(N ? N : 1));
        NAMELEN_CP = xmalloc(sizeof(size_t) * (size_t)(N ? N : 1));
        for (int s = 0; s < N; s++) {
            Str fsrc = fields[9 + s];
            char *d = xmalloc(fsrc.n + 1);
            size_t dn = 0;
            for (size_t i = 0; i < fsrc.n;) {
                if (i + 1 < fsrc.n && fsrc.p[i] == '.' && fsrc.p[i + 1] == '/') { i += 2; continue; }
                d[dn++] = fsrc.p[i++];
            }
            d[dn] = 0;
            NAMES[s].p = d; NAMES[s].n = dn;
            NAMELEN_CP[s] = cp_count((const unsigned char *)d, dn);
        }
        free(fields);
    }
    prof("header");
    if (N == 0) {
        buf_puts(&OUT, "\nSample names not found in VCF, your file may be corrupt or missing the header.\n\n");
        out_flush();
        return 0;
    }
    buf_puts(&OUT, "\nConverting file '"); buf_puts(&OUT, fn); buf_puts(&OUT, "':\n\n");
    buf_puts(&OUT, "Number of samples in VCF: "); buf_putll(&OUT, N); buf_putc(&OUT, '\n');

    if (OPT.msl > N) OPT.msl = N;
    /* output prefix */
    Buf pre = {0};
    if (OPT.prefix && OPT.prefix[0]) buf_puts(&pre, OPT.prefix);
    else {
        const char *base = strrchr(fn, '/');
        base = base ? base + 1 : fn;
        /* Path.name: last component (trailing slashes would make it a directory, not a file) */
        const char *p = base;
        int first = 1;
        for (;;) {
            const char *dot = strchr(p, '.');
            size_t len = dot ? (size_t)(dot - p) : strlen(p);
            if (len == 3 && tolower((unsigned char)p[0]) == 'v' && tolower((unsigned char)p[1]) == 'c' && tolower((unsigned char)p[2]) == 'f') break;
            if (!first) buf_putc(&pre, '.');
            buf_put(&pre, p, len);
            first = 0;
            if (!dot) break;
            p = dot + 1;
        }
    }
    buf_puts(&pre, ".min");
    buf_putll(&pre, OPT.msl);
    buf_putc(&pre, 0);
    if (stat(OPT.folder[0] ? OPT.folder : ".", &stt) != 0) mkdir_parents(OPT.folder);
    OUTFILE = path_join_norm(OPT.folder, pre.p);

    NUC = OPT.fasta || OPT.nexus || !OPT.phylipdisable;
    /* the reference creates its temporary files here; we create them only if it would leave them */
    Buf pth = {0};
    if (NUC) {   /* check that the temporary would be creatable, as open(...,"w") would */
        outfile_path(&pth, ".tmp");
        int fd = open(pth.p, O_WRONLY | O_CREAT | O_EXCL, 0644);
        if (fd < 0 && errno != EEXIST) { out_flush(); fprintf(stderr, "Traceback (most recent call last):\n  ...\nFileNotFoundError: [Errno %d] %s: '%s'\n", errno, strerror(errno), pth.p); exit(1); }
        if (fd >= 0) { close(fd); unlink(pth.p); }
    }
    if (OPT.write_used) {
        outfile_path(&pth, ".used_sites.tsv");
        fw_open(&USED, pth.p);
        buf_puts(&USED.b, "#CHROM\tPOS\tNUM_SAMPLES\n");
    }
    /* tile memory before spilling: SWIFTPHYLIP_MEM_MB, else 60% of MemAvailable (at least 256 MB) */
    size_t budget = 2048;
    {
        FILE *mi = fopen("/proc/meminfo", "r");
        if (mi) {
            char ln[256]; unsigned long long kb;
            while (fgets(ln, sizeof ln, mi)) if (sscanf(ln, "MemAvailable: %llu kB", &kb) == 1) { budget = (size_t)(kb * 6 / 10 / 1024); break; }
            fclose(mi);
        }
        if (budget < 256) budget = 256;
    }
    const char *mb = getenv("SWIFTPHYLIP_MEM_MB");
    if (mb && atol(mb) > 0) budget = (size_t)atol(mb);
    budget <<= 20;
    char *spill_dir = strdup(OUTFILE);
    { char *sl = strrchr(spill_dir, '/'); if (sl) { if (sl == spill_dir) sl[1] = 0; else *sl = 0; } else strcpy(spill_dir, "."); }
    if (NUC) mat_init(&NM, N, OPT.nexusbin ? budget / 2 : budget, spill_dir);
    if (OPT.nexusbin) mat_init(&BM, N, NUC ? budget / 2 : budget, spill_dir);

    /* pass 2: all lines */
    MAINC.is_main = 1; MAINC.nm = &NM; MAINC.bm = &BM;
    C = &MAINC;
    if (!gz && THREADS >= 2) parse_plain_parallel(fn, budget);
    else {
        Lines L; lines_open(&L, fn, gz, 256u << 10, 0);
        char *ls, *le;
        while (lines_next(&L, &ls, &le)) process_line(ls, le);
        prof("parse");
        lines_close(&L);
    }
    buf_puts(&OUT, "Total of genotypes processed: "); buf_putll(&OUT, MAINC.snp_num); buf_putc(&OUT, '\n');
    buf_puts(&OUT, "Genotypes excluded because they exceeded the amount of missing data allowed: "); buf_putll(&OUT, MAINC.snp_shallow); buf_putc(&OUT, '\n');
    buf_puts(&OUT, "Genotypes that passed missing data filter but were excluded for being MNPs: "); buf_putll(&OUT, MAINC.mnp_num); buf_putc(&OUT, '\n');
    buf_puts(&OUT, "SNPs that passed the filters: "); buf_putll(&OUT, MAINC.snp_accepted); buf_putc(&OUT, '\n');
    if (OPT.nexusbin) { buf_puts(&OUT, "Biallelic SNPs selected for binary NEXUS: "); buf_putll(&OUT, MAINC.snp_biallelic); buf_putc(&OUT, '\n'); }
    if (OPT.write_used) {
        buf_puts(&OUT, "Used sites saved to: '"); buf_puts(&OUT, OUTFILE); buf_puts(&OUT, ".used_sites.tsv'\n");
        fw_close(&USED);
    }
    buf_putc(&OUT, '\n');
    if (NUC) mat_finish(&NM);
    if (OPT.nexusbin) mat_finish(&BM);

    size_t longest = 0;
    for (int s = 0; s < N; s++) if (NAMELEN_CP[s] > longest) longest = NAMELEN_CP[s];
    int idx_out = -1;
    size_t ogl = strlen(og);
    for (int s = 0; s < N; s++) if (NAMES[s].n == ogl && memcmp(NAMES[s].p, og, ogl) == 0) { idx_out = s; break; }

    prof("finish");
    write_matrices(idx_out, longest);
    prof("write");

    /* stdout messages of the output phase, in the reference's order */
    if (idx_out >= 0) {
        if (NUC) { buf_puts(&OUT, "Outgroup, '"); buf_puts(&OUT, og); buf_puts(&OUT, "', added to the matrix(ces).\n"); }
        if (OPT.nexusbin) { buf_puts(&OUT, "Outgroup, '"); buf_puts(&OUT, og); buf_puts(&OUT, "', added to the binary matrix.\n"); }
    }
    for (int s = 0; s < N; s++) {
        if (s == idx_out) continue;
        if (NUC) {
            buf_puts(&OUT, "Sample "); buf_putll(&OUT, s + 1); buf_puts(&OUT, " of "); buf_putll(&OUT, N); buf_puts(&OUT, ", '");
            buf_put(&OUT, NAMES[s].p, NAMES[s].n); buf_puts(&OUT, "', added to the nucleotide matrix(ces).\n");
        }
        if (OPT.nexusbin) {
            buf_puts(&OUT, "Sample "); buf_putll(&OUT, s + 1); buf_puts(&OUT, " of "); buf_putll(&OUT, N); buf_puts(&OUT, ", '");
            buf_put(&OUT, NAMES[s].p, NAMES[s].n); buf_puts(&OUT, "', added to the binary matrix.\n");
        }
        out_maybe_flush();
    }
    buf_putc(&OUT, '\n');
    if (!OPT.phylipdisable) { buf_puts(&OUT, "PHYLIP matrix saved to: "); buf_puts(&OUT, OUTFILE); buf_puts(&OUT, ".phy\n"); }
    if (OPT.fasta) { buf_puts(&OUT, "FASTA matrix saved to: "); buf_puts(&OUT, OUTFILE); buf_puts(&OUT, ".fasta\n"); }
    if (OPT.nexus) { buf_puts(&OUT, "NEXUS matrix saved to: "); buf_puts(&OUT, OUTFILE); buf_puts(&OUT, ".nex\n"); }
    if (OPT.nexusbin) { buf_puts(&OUT, "BINARY NEXUS matrix saved to: "); buf_puts(&OUT, OUTFILE); buf_puts(&OUT, ".bin.nex\n"); }
    buf_puts(&OUT, "\nDone!\n\n");
    out_flush();
    return 0;
}
