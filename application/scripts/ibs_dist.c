/* Allele-sharing (IBS) distance from a vcf2phylip/SwiftPhylip binary NEXUS matrix (codes 0/1/2/?).
   D_ij = sum_{sites called in both} |g_i - g_j| / (2 * #sites called in both).
   Exact, via bit planes: |g_i-g_j| = [g_i>=1 != g_j>=1] + [g_i==2 != g_j==2].
   Usage: ibs_dist matrix.bin.nexus out_prefix   -> out_prefix.names.txt, out_prefix.dist.f64 (N*N doubles)
   Build: gcc -O3 -march=native -fopenmp ibs_dist.c -o ibs_dist */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <omp.h>

int main(int argc, char **argv) {
    if (argc < 3) { fprintf(stderr, "usage: ibs_dist matrix.bin.nexus out_prefix\n"); return 2; }
    FILE *f = fopen(argv[1], "rb");
    if (!f) { perror(argv[1]); return 1; }
    size_t cap = 1 << 20; char *line = malloc(cap); size_t len;
    long ntax = -1, nchar = -1; int inmat = 0;
    char **names = NULL; uint64_t *A = NULL, *B = NULL, *M = NULL; size_t W = 0; long n = 0; int anymiss = 0;
    ssize_t r;
    while ((r = getline(&line, &cap, f)) > 0) {
        len = (size_t)r; while (len && (line[len-1] == '\n' || line[len-1] == '\r')) line[--len] = 0;
        if (!inmat) {
            char *p = strstr(line, "NTAX="); if (p) ntax = atol(p + 5);
            p = strstr(line, "NCHAR="); if (p) nchar = atol(p + 6);
            if (strcmp(line, "MATRIX") == 0) {
                inmat = 1; W = (nchar + 63) / 64;
                names = calloc(ntax, sizeof(char*));
                A = calloc((size_t)ntax * W, 8); B = calloc((size_t)ntax * W, 8); M = calloc((size_t)ntax * W, 8);
                if (!A || !B || !M) { fprintf(stderr, "out of memory\n"); return 1; }
            }
            continue;
        }
        if (line[0] == ';') break;
        char *sp = strchr(line, ' '); char *last = strrchr(line, ' ');
        if (!sp) { fprintf(stderr, "bad matrix line\n"); return 1; }
        *sp = 0; names[n] = strdup(line);
        char *seq = last + 1; size_t L = strlen(seq);
        if ((long)L != nchar) { fprintf(stderr, "row %ld has %zu chars, expected %ld\n", n, L, nchar); return 1; }
        uint64_t *a = A + n * W, *b = B + n * W, *m = M + n * W;
        for (size_t k = 0; k < L; k++) {
            char c = seq[k]; uint64_t bit = 1ULL << (k & 63); size_t w = k >> 6;
            if (c == '?') { anymiss = 1; continue; }
            m[w] |= bit;
            if (c == '1' || c == '2') a[w] |= bit;
            if (c == '2') b[w] |= bit;
        }
        n++;
    }
    fclose(f);
    if (n != ntax) { fprintf(stderr, "read %ld rows, NTAX=%ld\n", n, ntax); return 1; }
    double *D = calloc((size_t)n * n, sizeof(double));
    const long BS = 32; const size_t CW = 512;
    long nb = (n + BS - 1) / BS;
    #pragma omp parallel for schedule(dynamic, 1)
    for (long bij = 0; bij < nb * nb; bij++) {
        long bi = bij / nb, bj = bij % nb; if (bj < bi) continue;
        uint64_t diff[32][32], cnt[32][32]; memset(diff, 0, sizeof diff); memset(cnt, 0, sizeof cnt);
        long i0 = bi * BS, j0 = bj * BS, i1 = i0 + BS < n ? i0 + BS : n, j1 = j0 + BS < n ? j0 + BS : n;
        for (size_t w0 = 0; w0 < W; w0 += CW) {
            size_t w1 = w0 + CW < W ? w0 + CW : W;
            for (long i = i0; i < i1; i++) {
                const uint64_t *ai = A + i * W, *bi_ = B + i * W, *mi = M + i * W;
                for (long j = (bi == bj ? i + 1 : j0); j < j1; j++) {
                    const uint64_t *aj = A + j * W, *bj_ = B + j * W, *mj = M + j * W;
                    uint64_t d = 0, c = 0;
                    if (!anymiss) {
                        for (size_t w = w0; w < w1; w++) d += __builtin_popcountll(ai[w] ^ aj[w]) + __builtin_popcountll(bi_[w] ^ bj_[w]);
                    } else {
                        for (size_t w = w0; w < w1; w++) { uint64_t mm = mi[w] & mj[w];
                            d += __builtin_popcountll((ai[w] ^ aj[w]) & mm) + __builtin_popcountll((bi_[w] ^ bj_[w]) & mm);
                            c += __builtin_popcountll(mm); }
                    }
                    diff[i - i0][j - j0] += d; cnt[i - i0][j - j0] += c;
                }
            }
        }
        for (long i = i0; i < i1; i++) for (long j = (bi == bj ? i + 1 : j0); j < j1; j++) {
            double c = anymiss ? (double)cnt[i - i0][j - j0] : (double)nchar;
            double v = c > 0 ? (double)diff[i - i0][j - j0] / (2.0 * c) : 0.0;
            D[i * n + j] = v; D[j * n + i] = v;
        }
    }
    char path[4096];
    snprintf(path, sizeof path, "%s.names.txt", argv[2]); FILE *o = fopen(path, "w");
    for (long i = 0; i < n; i++) fprintf(o, "%s\n", names[i]); fclose(o);
    snprintf(path, sizeof path, "%s.dist.f64", argv[2]); o = fopen(path, "wb");
    fwrite(D, sizeof(double), (size_t)n * n, o); fclose(o);
    fprintf(stderr, "ibs_dist: %ld samples, %ld sites, missing=%d\n", n, nchar, anymiss);
    return 0;
}
