# scfindME - adapting scfind to alternative splicing analysis

## Introduction
The emergence of large single cell RNA-seq(scRNA-seq) datasets facilitates biological discoveries at single-cell resolution. For single cell atlas to be useful to users from all disciplines, it should be organized in an accessible way. `scfind` builds an index for scRNA-seq datasets and allows fast searching with feature lists, and return cell types or individual cells [scfind](https://github.com/hemberg-lab/scfind).  
`scfindME` is built on top of `scfind` and it adaps `scfind` to alternative splicing analysis, thus revealing the co-regulated alternative splicing patterns across different cell types or individual cells.  
Quantification of single cell alternative splicing events are performed by [Whippet](https://github.com/timbitz/Whippet.jl) integreted in the workflow [MicroExonator](https://github.com/hemberg-lab/MicroExonator).

## Installation
To install or run `scfindME`, use the following codes in an R session:
```
install.packages("devtools")
devtools::install_github("YY-SONG0718/scfindME")
library(scfindME)
```
## User guide
Please refer to the `scfindME` package Vignettes for a detailed user guide.

## Contact
Please contact Miss Yuyao Song (ys6@sanger.ac.uk) for all ideas, questions, requests or bug reports.
