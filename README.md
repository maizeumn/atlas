This repository hosts code (primarily R and Python scripts) and data related to the maize developmental atlas dataset.

Repository architecture:
- `snk/`: pipeline scripts (snakemake and python)
- `src/`: data processing, statistical testing, visualization (R)
- `Rmd/`: figure and table typesetting (R markdown)
- `README.md`: (this file)

Demo datasets and scripts to:
- [obtain allele specific read counts](https://github.com/orionzhou/demo/tree/master/ase)
- [characterize regulatory patterns (cis/trans) using alelle-specific read counts](https://github.com/orionzhou/demo/blob/master/ase/cis_trans.md)

Here are links to a list of useful R scripts:
- [normalize raw read counts using TMM approach](/src/br.03.collect.R)
- [generate sample QC statistics, PCA, hierarchical clustering, etc.](/src/br.11.qc.sample.r)
- [run DE tests using DESeq2 or edgeR](/src/br.15.de.1.run.r)



