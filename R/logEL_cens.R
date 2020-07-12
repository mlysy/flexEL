#' Calculates log empirical likelihood for right-censored observations with optional continuity correction.
#'
#' @template args-omega
#' @template args-eps
#' @template args-delta
#' @template args-sp
#' @template args-support
#' @example examples/logEL_cens.R
#' @return Log empirical likelihood with continuity correction.
#' @export logEL_cens
logEL_cens <- function(omega, eps, delta, sp = 10, support = FALSE) {
  
  # input check
  if (sp <= 0) stop("s must be positive.")
  if (length(omega) != length(eps)) stop("omega must have the same length as eps.")
  
  if (sp == 0) {
    return(.LogELCens(omega, eps, delta, support))
  }
  else {
    return(.LogELSmooth(omega, eps, delta, sp, support))
  }
}
