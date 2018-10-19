require(tidyverse) 
require(rmarkdown) 
require(knitr) 
require(kableExtra) 
dirm = '~/projects/briggs/Rmd'
getwd()
rmarkdown::render(file.path(dirm, "rnaseq.r1.Rmd"))
#rmarkdown::render(file.path(dirm, "rnaseq.r1.m.Rmd"))

#rmarkdown::render(file.path(dirm, "rnaseq.figure.Rmd"))
#rmarkdown::render(file.path(dirm, "rnaseq.sm.Rmd"))
#rmarkdown::render(file.path(dirm, "rnaseq.unused.Rmd"))
