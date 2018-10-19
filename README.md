This repository hosts code (primarily R and Python scripts) and data related to the maize developmental atlas dataset.

Repository architecture:
- `snk/`: pipeline scripts (snakemake and python)
- `src/`: data processing, statistical testing, visualization (R)
- `Rmd/`: figure and table typesetting (R markdown)
- `README.md`: (this file)

In parcular, here are links to a list of R scripts potentially useful to the community:
- [normalize raw read counts using TMM approach](/src/br.03.collect.r)
- [generate sample QC statistics, PCA, hierarchical clustering, etc.](/src/br.11.qc.sample.r)
- [run DE tests using DESeq2 or edgeR](/src/br.15.de.1.run.r)
- [characterize additivity/dominance inheritance patterns](/src/br.15.de.2.r)
- [characterize regulatory patterns using alelle-specific read counts](/src/br.17.ase.1.run.r)



