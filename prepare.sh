#!/bin/bash
# Downloads the public inputs, builds the benchmark sets and builds SwiftPhylip.
# Needs Linux, gcc, make, zlib headers, curl and python3. About 300 MB of downloads and 200 MB of sets.
set -eu
cd "$(dirname "$0")"
mkdir -p data/raw data/dev data/heldout data/adversarial reference
B=http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502
VCF=ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
PANEL=integrated_call_samples_v3.20130502.ALL.panel
[ -s data/raw/$VCF ] || curl -sS --retry 10 -C - -o data/raw/$VCF $B/$VCF
[ -s data/raw/$PANEL ] || curl -sS --retry 5 -o data/raw/$PANEL $B/$PANEL
[ -s data/raw/pinf_sc50.vcf.gz ] || curl -sS --retry 5 -L -o data/raw/pinf_sc50.vcf.gz \
  https://raw.githubusercontent.com/knausb/pinfsc50/59d7324edda6d5875cf41797dc94a815ab2ad857/inst/extdata/pinf_sc50.vcf.gz
[ -s reference/vcf2phylip.py ] || curl -sS --retry 5 -L -o reference/vcf2phylip.py \
  https://raw.githubusercontent.com/edgardomortiz/vcf2phylip/85a260e655b5f3c5dc2cb9915eb24552a002357d/vcf2phylip.py
python3 -I scripts/prep_data.py pinf data/raw/pinf_sc50.vcf.gz data
(cd data && zcat raw/$VCF | python3 -I ../scripts/prep_data.py kg .)
python3 -I gen_data.py adversarial data/adversarial
sha256sum -c --quiet data/inputs.sha256
make -C solver
echo "prepared: solver/build/swiftphylip"
