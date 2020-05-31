#' Calculates log empirical likelihood given fully observed or right-censored observations.
#'
#' @template args-omega
#' @template args-eps
#' @template args-delta
#' @example examples/logEL.R
#' @return Log empirical likelihood
#' @export logEL
logEL <- function(omega, eps, delta, support = FALSE) {
  if (missing(eps) && missing(delta)) {
    logel <- .LogEL(omega, support)
  }
  else {
    logel <- .LogELC(omega,eps,delta,support)
  }
  return(logel)
}
