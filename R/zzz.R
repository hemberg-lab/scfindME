#' Essential File so module is loaded
#' 
#' @importFrom Rcpp loadModule
#' @useDynLib scfindME
#'
#' @include zzz.R
#' 
#' @export
NULL

loadModule('EliasFanoDB', TRUE)
