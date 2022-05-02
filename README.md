# scfindME - mining alternative splicing patterns in single-cell atlas

## Introduction
The emergence of large single cell RNA-seq(scRNA-seq) datasets facilitates biological discoveries at single-cell resolution. The information encoded in single cell atlases should be readily accessible to users from all disciplines. The Hemberg lab has developed `scfind`, a search engine for gene expression patterns in cell atlases [scfind](https://github.com/hemberg-lab/scfind).  

`scfindME` is built on top of `scfind` and it adapts `scfind` to **alternative splicing** analysis. Through querying with splicing features or cell types, `scfindME` conducts rapid searching and returns feature-enriched cell types or cell type-specific splicing signatures, respectively. `scfindME` also discovers cell type-specific mutually exclusive exon pairs and clustered spliced-in or spliced-out of a block of exons by enumerating over all possible combinations of splicing events. 

Quantification of single cell alternative splicing events is performed by [Whippet](https://github.com/timbitz/Whippet.jl) integrated in the workflow [MicroExonator](https://github.com/hemberg-lab/MicroExonator). scfindME builds custom indices for any single-cell dataset using the output files of `Whippet`.

## Installation
To install or run `scfindME`, use the following codes in an R session:
```
install.packages("devtools")
devtools::install_github("hemberg-lab/scfindME")
library(scfindME)
```
## User guide
Please refer to the `scfindME` package [Vignette](https://github.com/hemberg-lab/scfindME/blob/master/vignettes/example_workflow.ipynb) for a detailed user guide.

## Contact
Please contact Miss Yuyao Song (ys585@cam.ac.uk) for any enquiries or bug reports.
