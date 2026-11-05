#!/bin/bash
# Compares SwiftPhylip with vcf2phylip byte for byte, times both on the benchmark sets, and repeats the chromosome 22 analysis.
# Timings need about 10 minutes. The analysis needs python3 with numpy and matplotlib, gcc with OpenMP, about 5 GB of RAM and 6 GB of disk.
set -euo pipefail
cd "$(dirname "$0")"
python3 -I evaluate.py gate  --cand solver/build/swiftphylip
python3 -I evaluate.py final --cand solver/build/swiftphylip
[ "${SKIP_APPLICATION:-0}" = "1" ] && exit 0
ROOT=$PWD
VCF=ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
mkdir -p application/out
(cd application/out && "$ROOT/solver/build/swiftphylip" -i "$ROOT/data/raw/$VCF" -b --output-folder conv > stdout.txt)
gcc -O3 -march=native -fopenmp application/scripts/ibs_dist.c -o application/scripts/ibs_dist
application/scripts/ibs_dist "application/out/conv/${VCF%.vcf.gz}.min4.bin.nexus" application/out/d
python3 -I application/scripts/nj_analysis.py application/out/d data/raw/integrated_call_samples_v3.20130502.ALL.panel application/out
python3 -I application/scripts/explore_clades.py application/out
