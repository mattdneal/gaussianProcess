# The first two lines of the roxygen code below create the entries in NAMESPACE
# required by Rcpp.

#' Utilities for training Gaussian processes.
#'
#' gaussianProcess provides functionality to define kernels (either
#' stand-alone or as a combination of basic kernels) and also
#' to automatically perform model selection from a basic set of kernels
#' to select the optimum kernel for modelling the input data, using the
#' Bayesian Information Criterion as a search criteria.
#'
#'@importFrom Rcpp evalCpp
#'@useDynLib gaussianProcess
#'
"_PACKAGE"
