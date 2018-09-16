#' Calculate partial sum of omegas with smoothing
#'
#' @param omegas Vector of probability weights
#' @param epsilons Vector of residuals.
#' @param ii An integer indicator between 1 and the length of omegas/epsilons.
#' @return Partial sum of omegas according to residuals that are no larger than the ii-th residual.
#' @details ...
#' @export evalPsos.smooth
evalPsos.smooth <- function(ii, omegas, epsilons, s=10) {
  # input check here
  if (s <=0) stop("s must be positive.")
  if (length(omegas) != length(epsilons)) stop("omegas must have the same length as epsilons.")
  if (ii > length(omegas) || ii <= 0) stop("index must be 1 to nObs.")
  return(.evalPsos.Smooth(ii-1, omegas, epsilons, s))
}
