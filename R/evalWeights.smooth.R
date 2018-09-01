#' Calculate partial sum of omegas with smoothing
#'
#' @param deltas Vector of censroing indicators.
#' @param omegas Vector of probability weights
#' @param epsilons Vector of residuals.
#' @return Partial sum of omegas according to residuals that are no larger than the ii-th residual.
#' @details ...
#' @export evalPsos.smooth
evalWeights.smooth <- function(deltas, omegas, epsilons, s=10) {
  # input check here
  if (s <=0) stop("s must be positive.")
  if (length(omegas) != length(epsilons)) stop("omegas must have the same length as epsilons.")
  if (length(omegas) != length(deltas)) stop("omegas must have the same length as deltas.")
  return(.evalWeights.smooth(deltas, omegas, epsilons, s))
}
