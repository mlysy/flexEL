#' Calculate the smoothed log empirical likelihood under censoring
#'
#' @template args-omega
#' @template args-eps
#' @template args-sp
#' @return Smoothed log empirical likelihood (a scalar).
#' @details ...
#' @export logEL.smooth
logEL.smooth <- function(omega, eps, deltas, sp=10, support = FALSE) {
  # input check here
  if (sp <=0) stop("s must be positive.")
  if (length(omega) != length(eps)) stop("omega must have the same length as eps.")
  return(.logEL.smooth(omega,eps,deltas,sp,support))
}
