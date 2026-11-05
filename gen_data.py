#!/usr/bin/env python3
"""Adversarial and random VCF generator for the SwiftPhylip correctness gate.

gen_case(rng) returns one random case: (input file name, input bytes, argument list).
Cases cover every output-affecting vcf2phylip option (-p -f -n -b -m -o -w --output-prefix
--output-folder), plain / gzip / multi-member gzip input, odd file names, any ploidy, missing and
half-missing calls, multi-allelic and odd ALT fields, GATK <NON_REF>, deletions (* and -),
lowercase and IUPAC alleles, malformed and anomalous lines, CRLF and lone-CR line ends, blank,
comment and whitespace-padded lines, Unicode whitespace at line ends, UTF-8 sample names, and the
inputs that make vcf2phylip raise (exit status 1 with partial temporary files).

Input domain: the decompressed input is valid UTF-8 text. (vcf2phylip's behaviour on invalid
UTF-8 depends on its internal decode-buffer boundaries.) The -r option (random resolution) and
-h / -v are outside the target.

`python3 gen_data.py adversarial OUTDIR` writes the fixed adversarial files used by the gate.
"""
import gzip
import io
import random
import sys
import os

NAME_CHARS = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789_-."
UTF8_CHARS = ["é", "ß", "中", "😀", "ø", "Å", "ı", " x"]
UWS = ["\x0b", "\x0c", "\x1c", "\x1d", "\x1e", "\x1f", "\x85", "\xa0", "　", " ", " ", "\t"]


def rand_name(rng):
    r = rng.random()
    if r < 0.03:
        return ""
    n = rng.randint(1, 14) if rng.random() < 0.9 else rng.randint(15, 40)
    s = "".join(rng.choice(NAME_CHARS) for _ in range(n))
    if rng.random() < 0.08:
        p = rng.randint(0, len(s))
        s = s[:p] + "./" + s[p:]
    if rng.random() < 0.08:
        p = rng.randint(0, len(s))
        s = s[:p] + rng.choice(UTF8_CHARS) + s[p:]
    if rng.random() < 0.03:
        s = s + " " + s[:2]
    return s


CLEAN_REF = ["A", "C", "G", "T", "a", "c", "g", "t", "N", "n", "-", "*", "AT", "ACG", "<NON_REF>"]
CLEAN_ALT = {"A", "C", "G", "T", "A,C", "C,G,T", "A,T", "*", "-", "A,*", ".", "<NON_REF>", "A,<NON_REF>",
             "<NON_REF>,C", "AT", "A,TT", "<DEL>", "<CN0>", "N", "A,N", "A,C,G,T,A,C,G,T,A,C,G", "-,A", "T,G,C,A,*"}


def rand_ref(rng, wild=True):
    if not wild:
        return rng.choice(CLEAN_REF[:4] * 6 + CLEAN_REF)
    return rng.choices(
        ["A", "C", "G", "T", "a", "c", "g", "t", "N", "n", "-", "*", "AT", "ACG", "", "R", "M", ".", "é", "<NON_REF>"],
        [15, 15, 15, 15, 2, 2, 2, 2, 2, 1, 2, 2, 6, 2, 1, 1, 1, 1, 0.5, 0.5])[0]


def rand_alt(rng, wild=True):
    if not wild:
        while True:
            a = rand_alt(rng, True)
            if a in CLEAN_ALT:
                return a
    return rng.choices(
        ["A", "C", "G", "T", "A,C", "C,G,T", "A,T", "g", "t,A", "*", "-", "A,*", ".", "<NON_REF>",
         "A,<NON_REF>", "<NON_REF>,C", "AT", "A,TT", "<DEL>", "<CN0>", ",AC", "AC,", "A,,CG", "N", "A,N",
         "A,C,G,T,A,C,G,T,A,C,G", "", "-,A", "T,G,C,A,*"],
        [12, 12, 12, 12, 5, 3, 3, 1.5, 1, 2, 2, 2, 3, 3, 3, 1, 4, 1.5, 1, 1, 0.7, 0.7, 0.7, 1.5, 1, 0.8, 0.5, 1, 1])[0]


def rand_gt(rng, nalleles, ploidy, miss, wild=True):
    r = rng.random()
    if not wild and r < 0.012:
        r = 0.5
    if not wild and nalleles == 0:
        return rng.choice(["0/0", "./.", "0|0", ".", "0", "0/0:3"])
    if r < 0.012:
        return rng.choice(["", ":", "/", "|"])                   # empty genotype -> vcf2phylip raises
    if r < 0.03:
        return rng.choice(["a/b", "0/" + str(min(nalleles + 1, 9)), "x", "0/ 1", "1/-", "0/1 "])  # malformed
    if r < 0.05:
        return rng.choice([".", "./.", ".|.", "0/.", "./1", ".:5", "1/.:3", "0//1", "01", "10", "0|0|1"] + (["0/1/2/3/."] if wild else []))
    if r < 0.06 and nalleles > 1:
        return "%d/%d" % (rng.randint(0, nalleles - 1), rng.randint(0, nalleles - 1)) + "0"
    parts = []
    for _ in range(ploidy):
        if rng.random() < miss:
            parts.append(".")
        else:
            parts.append(str(rng.randint(0, max(0, nalleles - 1))))
    sep = rng.choice(["/", "|"]) if rng.random() < 0.95 else rng.choice(["/", "|"])
    g = sep.join(parts)
    if rng.random() < 0.3:
        g += ":" + str(rng.randint(0, 99)) + rng.choice(["", ":3,4", ":.:1"])
    return g


def gen_vcf(rng, wild):
    ns = rng.choices([0, 1, 2, 3, 4, 5, 6, 8, 10, 12, 20, 40, 60], [2, 4, 6, 8, 10, 10, 10, 10, 10, 8, 6, 3, 2])[0]
    names = [rand_name(rng) for _ in range(ns)]
    if ns > 1 and rng.random() < 0.08:
        names[rng.randrange(ns)] = names[rng.randrange(ns)]      # duplicate name
    lines = []
    lines.append("##fileformat=VCFv4.2")
    for _ in range(rng.randint(0, 3)):
        lines.append("##INFO=<ID=X%d,Number=1,Type=Integer,Description=\"d\">" % rng.randint(0, 9))
    hdr_mode = rng.random()
    header = "\t".join(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"] + names)
    if hdr_mode < 0.02:
        pass                                                     # no #CHROM line
    elif hdr_mode < 0.03:
        lines.append(" " + header)                               # not recognised as header
    else:
        lines.append(header)
        if rng.random() < 0.02:
            lines.append(header + "\textra")                     # second #CHROM line
    miss = rng.choice([0.0, 0.0, 0.05, 0.2, 0.5, 0.9])
    nrec = rng.choices([0, 1, 3, 8, 20, 40, 80], [2, 4, 8, 20, 30, 25, 11])[0]
    chrom = rng.choice(["1", "chr22", "Supercontig_1.50", "scaf_é", "X"])
    pos = rng.randint(1, 1000)
    for _ in range(nrec):
        pos += rng.randint(0, 50)
        ref = rand_ref(rng, wild)
        alt = rand_alt(rng, wild)
        nall = 1 + len(alt.split(","))
        ploidy = rng.choices([1, 2, 3, 4, 5], [8, 80, 5, 5, 2 if wild else 0])[0]
        gts = [rand_gt(rng, 0 if (alt == "." and not wild) else nall, ploidy, miss, wild) for _ in range(ns)]
        fields = [chrom, str(pos), rng.choice([".", "rs%d" % pos]), ref, alt,
                  rng.choice([".", "50", "3.2"]), rng.choice([".", "PASS"]),
                  rng.choice([".", "DP=10", "AC=1;AN=2"]), rng.choice(["GT", "GT:DP", "GT:AD:DP"])] + gts
        r = rng.random()
        if r < 0.03:
            fields = fields[:-1] if len(fields) > 9 else fields + ["0/0"]   # anomalous field count
        elif r < 0.04:
            fields = fields + ["0/1"]
        line = "\t".join(fields)
        r = rng.random()
        if r < 0.03:
            line = line + rng.choice([" ", "  ", "\t", " \t"])
        elif r < 0.05:
            line = rng.choice([" ", "\t", "  "]) + line
        elif r < 0.065:
            line = line + rng.choice(UWS)
        elif r < 0.075:
            line = rng.choice(UWS[:9]) + line
        lines.append(line)
        r = rng.random()
        if r < 0.04:
            lines.append(rng.choice(["", " ", "\t", "  \t "]))
        elif r < 0.06:
            lines.append(rng.choice(["#comment", "  #indented comment", "##late meta"]))
    eol = rng.choices(["\n", "\r\n", "mixed"], [85, 11, 4])[0]
    out = []
    for ln in lines:
        if eol == "mixed":
            e = rng.choice(["\n", "\r\n", "\r", "\n"])
            if rng.random() < 0.1 and len(ln) > 4:
                p = rng.randint(1, len(ln) - 1)
                ln = ln[:p] + "\r" + ln[p:]                      # lone CR splits a line
        else:
            e = eol
        out.append(ln + e)
    text = "".join(out)
    if rng.random() < 0.1:
        text = text.rstrip("\r\n")                               # no final newline
    return names, text.encode("utf-8")


def compress(rng, data):
    mode = rng.choices(["plain", "gz", "multi"], [60, 25, 15])[0]
    if mode == "plain":
        return data, False
    if mode == "gz":
        return gzip.compress(data, compresslevel=rng.randint(1, 9), mtime=0), True
    bio = io.BytesIO()
    pos = 0
    while pos < len(data) or pos == 0:
        step = rng.randint(1, max(1, len(data) // 3 + 1))
        bio.write(gzip.compress(data[pos:pos + step], mtime=0))
        pos += step
        if pos >= len(data):
            break
    return bio.getvalue(), True


def gen_case(rng):
    wild = rng.random() < 0.25          # 25% of cases allow inputs that make vcf2phylip raise
    names, text = gen_vcf(rng, wild)
    data, gz = compress(rng, text)
    base = rng.choice(["sample", "my.data", "x", ".hidden", "a.vcf.b", "data", "UPPER.set", "run_1"])
    ext = rng.choice([".vcf", ".vcf", ".VCF", ""])
    if gz:
        ext += rng.choice([".gz", ".gz", ".GZ"])
    fname = base + ext
    args = []
    for flag, p in (("-p", 0.2), ("-f", 0.4), ("-n", 0.4), ("-b", 0.4), ("-w", 0.3)):
        if rng.random() < p:
            args.append(rng.choice([flag, {"-p": "--phylip-disable", "-f": "--fasta", "-n": "--nexus",
                                            "-b": "--nexus-binary", "-w": "--write-used-sites"}[flag]]))
    if rng.random() < 0.15 and len(args) >= 2 and all(a.startswith("-") and len(a) == 2 for a in args[:2]):
        args = ["-" + args[0][1] + args[1][1]] + args[2:]       # combined short flags
    if rng.random() < 0.5:
        m = rng.choice([-2, 0, 1, 2, 3, 4, 5, len(names), len(names) + 3, 10 ** 6])
        args += rng.choice([["-m", str(m)], ["-m%d" % m], ["--min-samples-locus=%d" % m],
                            ["--min-samples-locus", str(m)]])
    r = rng.random()
    if r < 0.5:
        pass
    else:
        if r < 0.75 and names:
            og = rng.choice(names)
        elif r < 0.85:
            og = "not_a_sample"
        elif r < 0.9 and names:
            og = rng.choice(names) + "," + rng.choice(names)
        elif r < 0.95 and names:
            og = rng.choice(names) + ";zz"
        else:
            og = ""
        args += ["--outgroup=" + og] if (og.startswith("-") or og == "" or rng.random() < 0.3) else ["-o", og]
    if rng.random() < 0.3:
        args += ["--output-prefix", rng.choice(["pre", "p.q", "with space", "Z", "out.vcf"])]
    r = rng.random()
    if r < 0.5:
        pass
    else:
        args += ["--output-folder", rng.choice(["out", "out/sub", "./o2/", "o3//x", "."])]
    if rng.random() < 0.05 and args:
        args = args + [args[0]] if args[0] in ("-f", "-n", "-b", "-w", "-p") else args
    return fname, data, args


def adversarial(outdir):
    os.makedirs(outdir, exist_ok=True)
    rng = random.Random(20261104)
    files = {}
    # 1. > 500000 records to exercise the progress message, 2 samples, with anomalous and blank lines.
    buf = ["##fileformat=VCFv4.2", "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2"]
    for i in range(520000):
        if i % 50000 == 7:
            buf.append("1\t%d\t.\tA\tC\t.\t.\t.\tGT\t0/1" % i)            # anomalous (one sample)
        elif i % 70001 == 3:
            buf.append("")
        else:
            buf.append("1\t%d\t.\t%s\t%s\t.\t.\t.\tGT\t%s\t%s" % (
                i, "ACGT"[i % 4], "CGTA"[i % 4] if i % 9 else "G,T", "0/1" if i % 3 else "1|1", "./." if i % 11 == 0 else "0/0"))
    files["adv_progress.vcf"] = ("\n".join(buf) + "\n").encode()
    # 2. wide: 3000 samples with long and UTF-8 names, mixed ploidy, 30 records, CRLF, gzip multi-member.
    names = ["smp%04d" % i + ("é" * (i % 5)) + ("中" if i % 97 == 0 else "") for i in range(3000)]
    buf = ["##fileformat=VCFv4.2", "\t".join(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"] + names)]
    for i in range(30):
        alt = ["C", "C,G", "*", "<NON_REF>", "A,C,T"][i % 5]
        nall = 1 + len(alt.split(","))
        gts = [rand_gt(rng, nall, [2, 2, 1, 3][j % 4], 0.1) for j in range(3000)]
        gts = [g if g not in ("", ":", "/", "|") else "./." for g in gts]
        buf.append("\t".join(["chr2", str(100 + i), ".", "ACGT"[i % 4], alt, ".", ".", ".", "GT"] + gts))
    data = ("\r\n".join(buf) + "\r\n").encode()
    bio = io.BytesIO()
    for k in range(0, len(data), 65280):
        bio.write(gzip.compress(data[k:k + 65280], mtime=0))
    files["adv_wide_crlf.vcf.gz"] = bio.getvalue()
    # 3. medium random file: 200 samples x 3000 records from the fuzz generator's record model.
    names = ["ind_%d" % i for i in range(200)]
    buf = ["##fileformat=VCFv4.2", "\t".join(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"] + names)]
    for i in range(3000):
        ref = rand_ref(rng)
        alt = rand_alt(rng)
        if alt in ("g", "t,A", "N", "A,N", "A,C,G,T,A,C,G,T,A,C,G", "T,G,C,A,*") or ref in ("é", "R", "M", "<NON_REF>", ".", "a", "c", "g", "t", "n"):
            alt = "A"
        nall = 1 + len(alt.split(","))
        gts = [rand_gt(rng, nall, 2, 0.3) for _ in range(200)]
        gts = [g if g not in ("", ":", "/", "|", "0/1/2/3/.") else "./." for g in gts]
        buf.append("\t".join(["chr5", str(1000 + i * 7), ".", ref, alt, ".", ".", ".", "GT"] + gts))
    files["adv_medium.vcf"] = ("\n".join(buf) + "\n").encode()
    for k, v in files.items():
        with open(os.path.join(outdir, k), "wb") as f:
            f.write(v)
        print(k, len(v))


if __name__ == "__main__":
    if len(sys.argv) >= 3 and sys.argv[1] == "adversarial":
        adversarial(sys.argv[2])
    else:
        rng = random.Random(int(sys.argv[1]) if len(sys.argv) > 1 else 1)
        f, d, a = gen_case(rng)
        print(f, len(d), a)
