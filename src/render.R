#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("argparse"))

parser <- ArgumentParser(description = 'render .Rmd file')
parser$add_argument("fi", help = "input Rmd file")
args <- parser$parse_args()
fi = args$fi
stopifnot(file.exists(fi))

require(tidyverse) 
require(rmarkdown) 
require(knitr) 
require(kableExtra) 

rmarkdown::render(fi, clean = T)
