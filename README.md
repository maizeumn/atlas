This repository hosts code (primarily R and Python scripts) and data related to the maize developmental atlas dataset.

Repository architecture:
- `snk/`: pipeline scripts (snakemake and python)
- `src/`: data processing, statistical testing, visualization (R)
- `Rmd/`: figure and table typesetting (R markdown)
- `README.md`: (this file)

Demo datasets and scripts to:
- [obtain allele specific read counts](https://github.com/orionzhou/demo/tree/master/ase)
- [cis_trans.R](https://github.com/orionzhou/demo/blob/master/cis_trans/README.md): Classify cis/trans inheritance pattern using inbred/hybrid expression, in `basic` mode or `differential` mode
  ```
    $ ./cis_trans.R -h
    usage: ./cis_trans.R [-h] [--mode MODE] [--min_rc MIN_RC] [--n_cpu N_CPU]
                         f_rc f_sf f_dsp fo

    Classify cis/trans inheritance pattern using inbred/hybrid RNA-Seq read counts

    positional arguments:
      f_rc             read count table
      f_sf             sample-wise size factor table
      f_dsp            gene-wise dispersion table
      fo               output file

    optional arguments:
      -h, --help       show this help message and exit
      --mode MODE      cis/trans test mode, "basic" for steady-state cis/trans
                       test and "diff" for control/treatment differential test
                       [default: basic]
      --min_rc MIN_RC  minimum read counts to filter low-expressed genes [default:
                       10]
      --n_cpu N_CPU    number of CPUs / threads to use for parallel processing
                       (for spped up if you have many genes) [default: 1]
  ```

Here are links to a list of useful R scripts:
- [normalize raw read counts using TMM approach](/src/br.03.collect.R)
- [generate sample QC statistics, PCA, hierarchical clustering, etc.](/src/br.11.qc.sample.r)
- [run DE tests using DESeq2 or edgeR](/src/br.15.de.1.run.r)



