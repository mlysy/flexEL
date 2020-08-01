#' Fast and flexible empirical likelihood calculations.
#' 
#' A flexible framework for Empirical Likelihood methods for regression problems. The bulk of computations are done in C++ for speed considerations.  The C++ code is exposed as a header-only library.  The package provides mean and quantile regression with (optionally) right-censored responses using location and location-scale models, but users can added in any regression models as C++ header files or use the R interface directly to solve their problems.
#'
#' @docType package
#' @name flexEL
#' @useDynLib flexEL, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @importFrom Rdpack reprompt
#' @keywords internal
NULL
