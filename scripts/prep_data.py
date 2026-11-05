#!/usr/bin/env python3
"""Slice public real VCFs into development and held-out sets (chosen before any optimisation).
Usage: zcat chr22.vcf.gz | prep_data.py kg OUTDIR ; prep_data.py pinf pinf_sc50.vcf.gz OUTDIR"""
import sys, gzip, os

def write(path, header_meta, chrom_fields, sample_idx, records):
    gz = path.endswith(".gz")
    f = gzip.open(path, "wb", compresslevel=6) if gz else open(path, "wb")
    for h in header_meta: f.write(h)
    cols = chrom_fields[:9] + [chrom_fields[9 + i] for i in sample_idx]
    f.write(b"\t".join(cols) + b"\n")
    for r in records:
        fs = r.rstrip(b"\n").split(b"\t")
        f.write(b"\t".join(fs[:9] + [fs[9 + i] for i in sample_idx]) + b"\n")
    f.close()

def kg(outdir):
    # (name, first record index, n records, sample selector)
    specs = [
        ("dev/kg_all2504_r100k.vcf",      100000, 4000,  lambda n: list(range(n))),
        ("dev/kg_s250_r300k.vcf.gz",      300000, 40000, lambda n: list(range(0, n, 10))[:250]),
        ("heldout/kg_all2504_r700k.vcf",  700000, 4000,  lambda n: list(range(n))),
        ("heldout/kg_s500_r900k.vcf.gz",  900000, 20000, lambda n: list(range(2, n, 5))[:500]),
    ]
    meta, chrom = [], None
    buf = {s[0]: [] for s in specs}
    idx = -1
    last = max(s[1] + s[2] for s in specs)
    for line in sys.stdin.buffer:
        if line.startswith(b"##"): meta.append(line); continue
        if line.startswith(b"#CHROM"): chrom = line.rstrip(b"\n").split(b"\t"); continue
        idx += 1
        for name, a, n, _ in specs:
            if a <= idx < a + n: buf[name].append(line)
        if idx >= last: break
    ns = len(chrom) - 9
    for name, a, n, sel in specs:
        write(os.path.join(outdir, name), meta, chrom, sel(ns), buf[name])
        print(name, len(buf[name]), len(sel(ns)))

def pinf(src, outdir):
    meta, chrom, recs = [], None, []
    with gzip.open(src, "rb") as f:
        for line in f:
            if line.startswith(b"##"): meta.append(line)
            elif line.startswith(b"#CHROM"): chrom = line.rstrip(b"\n").split(b"\t")
            else: recs.append(line)
    half = len(recs) // 2
    ns = len(chrom) - 9
    write(os.path.join(outdir, "dev/pinf_sc50_first_half.vcf.gz"), meta, chrom, list(range(ns)), recs[:half])
    write(os.path.join(outdir, "heldout/pinf_sc50_second_half.vcf"), meta, chrom, list(range(ns)), recs[half:])
    print("pinf", half, len(recs) - half, ns)

if __name__ == "__main__":
    if sys.argv[1] == "kg": kg(sys.argv[2])
    else: pinf(sys.argv[2], sys.argv[3])
