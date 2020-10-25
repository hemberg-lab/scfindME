# scfindME - a search engine for alternative splicing patterns in single-cell atlas

## Introduction
The emergence of large single cell RNA-seq(scRNA-seq) datasets facilitates biological discoveries at single-cell resolution. For single cell atlas to be useful to users from all disciplines, it should be organized in an accessible way. `scfind` builds an index for scRNA-seq datasets and allows fast searching with feature lists, and return cell types or individual cells [scfind](https://github.com/hemberg-lab/scfind).  
`scfindME` is built on top of `scfind` and it adaps `scfind` to alternative splicing analysis. Through receiving queries of splicing features or cell types, `scfindME` conducts rapid searching and returns intuitive results. Thus, `scfindME` facilitates the discovery of co-regulated alternative splicing patterns in large-scale single-cel atlas.  
Quantification of single cell alternative splicing events are performed by [Whippet](https://github.com/timbitz/Whippet.jl) integreted in the workflow [MicroExonator](https://github.com/hemberg-lab/MicroExonator). scfindME builds custom indices for any single-cell dataset using the output files of `Whippet`.

## Installation
To install or run `scfindME`, use the following codes in an R session:
```
install.packages("devtools")
devtools::install_github("hemberg-lab/scfindME")
library(scfindME)
```
## User guide
Please refer to the `scfindME` package Vignettes for a detailed user guide.

## Contact
Please contact Miss Yuyao Song (ys6@sanger.ac.uk) for all ideas, questions, requests or bug reports.
