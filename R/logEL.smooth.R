#' Calculate the smoothed log empirical likelihood under censoring
#'
#' @param omegas Vector of probability weights
#' @param epsilons Vector of residuals.
#' @return A scalar of log empirical likelihood
#' @details ...
#' @export logEL.smooth
logEL.smooth <- function(omegas, epsilons, deltas, s=10) {
  # input check here
  if (s <=0) stop("s must be positive.")
  if (length(omegas) != length(epsilons)) stop("omegas must have the same length as epsilons.")
  return(.logEL.smooth(omegas,epsilons,deltas,s))
}
