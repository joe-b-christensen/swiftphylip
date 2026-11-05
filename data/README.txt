prepare.sh downloads the public inputs into data/raw/ and reference/ and builds the benchmark sets. data/inputs.sha256 holds
checksums for the downloads and the uncompressed sets; gzip files are not listed because their bytes depend on the zlib build.

Downloads:
- 1000 Genomes Phase 3 chromosome 22, v5b (205,612,353 bytes) and the sample panel:
  http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/
- pinfsc50 (Knaus and Grunwald 2017), commit 59d7324:
  https://github.com/knausb/pinfsc50
- vcf2phylip 2.9, commit 85a260e:
  https://github.com/edgardomortiz/vcf2phylip

Sets, fixed before development (scripts/prep_data.py; 0-based record indices of the chromosome 22 file):
  dev/kg_all2504_r100k.vcf        records 100000-103999, all 2504 samples
  dev/kg_s250_r300k.vcf.gz        records 300000-339999, samples 0,10,...,2490 (250)
  dev/pinf_sc50_first_half.vcf.gz pinfsc50 records 1-11015 (18 isolates)
  heldout/kg_all2504_r700k.vcf    records 700000-703999, all 2504 samples
  heldout/kg_s500_r900k.vcf.gz    records 900000-919999, samples 2,7,...,2497 (500)
  heldout/pinf_sc50_second_half.vcf pinfsc50 records 11016-22031
  adversarial/                    three fixed generated files (gen_data.py, seed 20261104)
