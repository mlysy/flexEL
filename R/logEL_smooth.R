#' Calculate log empirical likelihood under censoring with continuity correction
#'
#' @template args-omega
#' @template args-eps
#' @template args-sp
#' @example examples/logEL_smooth.R
#' @return Log empirical likelihood with continuity correction
#' @export logEL_smooth
logEL_smooth <- function(omega, eps, deltas, sp=10, support = FALSE) {
  # input check here
  if (sp <= 0) stop("s must be positive.")
  if (length(omega) != length(eps)) stop("omega must have the same length as eps.")
  return(.LogELSmooth(omega,eps,deltas,sp,support))
}
