# Application

A neighbour-joining tree of the 2,504 individuals of the 1000 Genomes Project from chromosome 22 SNPs, and how well its groups
agree with the five super-population labels. This covers one chromosome only.

`reproduce.sh` converts the VCF with SwiftPhylip (`-b`), computes allele-sharing distances (`scripts/ibs_dist.c`), builds the
tree and the agreement measures (`scripts/nj_analysis.py`) and looks at the clades outside the largest pure clade
(`scripts/explore_clades.py`; exploratory). `SKIP_APPLICATION=1 ./reproduce.sh` leaves this step out.

`results/` holds the published results: `application_results.json`, `exploratory_clades.json` and the tree, `tree.nwk`.
