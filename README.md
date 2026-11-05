# swiftphylip

swiftphylip is a drop-in replacement for [vcf2phylip](https://github.com/edgardomortiz/vcf2phylip) 2.9. It converts a VCF
into the same PHYLIP, FASTA, NEXUS and binary NEXUS files, prints the same messages and exits with the same status, for every
deterministic option.

On six real data sets on an AVX2 laptop, it was 12- to 706-fold faster than vcf2phylip. The gain grows with the number of
samples, because vcf2phylip re-reads its temporary file once per sample.

If you need only PHYLIP from SNPs, [plink2](https://www.cog-genomics.org/plink/2.0/) is also a good choice:
`plink2 --vcf in.vcf.gz --snps-only just-acgt --export phylip` wrote the same PHYLIP bytes as vcf2phylip on all four 1000 Genomes
sets tested, and on two *P. infestans* sets with `--geno 0.7778` added in place of vcf2phylip's `-m 4`. swiftphylip was 2.5- to
5-fold faster than plink2 on those six sets, and under 2-fold faster on a slice of 200,000 records from 2,504 samples. plink2 doesn't
write FASTA or NEXUS and doesn't accept every allele and ploidy that vcf2phylip does.

## Build

Requirements: Linux, gcc or clang, make and zlib headers.

```
make -C solver         # produces solver/build/swiftphylip
```

AVX2 and AVX-512 are used when the CPU has them, with a scalar fallback. All paths write the same bytes.

## Usage

```
solver/build/swiftphylip -i calls.vcf.gz                # PHYLIP
solver/build/swiftphylip -i calls.vcf.gz -f -n -b       # PHYLIP, FASTA, NEXUS and binary NEXUS
```

- Options work as in vcf2phylip: `-i`, `--output-folder`, `--output-prefix`, `-m`, `-o`, `-p`, `-f`, `-n`, `-b`, `-w`.
- `-r`, which resolves heterozygotes at random, is refused.
- `SWIFTPHYLIP_THREADS` sets the number of threads (default: all CPUs). `SWIFTPHYLIP_MEM_MB` sets the memory held before the matrix
  spills to a temporary file (default: 60% of available memory).
- Input that isn't valid UTF-8 was not tested.

## Check it yourself

```
./prepare.sh           # downloads the public data, builds the benchmark sets and swiftphylip
./reproduce.sh         # byte comparison, timings and the chromosome 22 analysis
```

`evaluate.py gate` compares every output file, leftover temporary file, standard output and exit status with vcf2phylip on fresh
random inputs, three fixed inputs and the development sets, and prints `GATE PASS` when all of them match. `reproduce.sh` also
needs python3 with numpy and matplotlib and gcc with OpenMP. Its last step, the chromosome 22 analysis, needs about 5 GB of RAM and
6 GB of disk; `SKIP_APPLICATION=1 ./reproduce.sh` leaves it out.

`data/README.txt` gives the sources of the data, and `application/` holds the analysis scripts and results.

## Licence

GPL-3.0 (`LICENSE`), as vcf2phylip. swiftphylip reimplements vcf2phylip's behaviour.

## Citation

See `CITATION.cff`. Please also cite vcf2phylip (Ortiz EM, Zenodo, 2019, doi:10.5281/zenodo.2540861).
