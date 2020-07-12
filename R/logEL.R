#' Calculates log empirical likelihood given fully observed observations.
#'
#' @template args-omega
#' @example examples/logEL.R
#' @return Log empirical likelihood.
#' @export logEL
logEL <- function(omega, support = FALSE) {
  return(.LogEL(omega, support))
}
