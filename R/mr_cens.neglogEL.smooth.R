#' Calculate the negative log EL with smoothing for mean regression under censoring
#'
#' @param omegas Vector of probability weights
#' @param epsilons Vector of residuals.
#' @param ii An integer indicator between 1 and the length of omegas/epsilons.
#' @return Partial sum of omegas according to residuals that are no larger than the ii-th residual.
#' @details ...
#' @export mr_cens.neglogEL.smooth
mr_cens.neglogEL.smooth <- function(y, X, deltas, beta, s=10) {
  G <- mr.evalG(y, X, beta)
  epsilons <- evalEpsilons(y,X,beta)
  omegas <- omega.hat.EM.smooth(G,deltas,epsilons,s)
  if (!anyNA(omegas)) {
    res <- -logEL.smooth(omegas,epsilons,deltas)
    # gradlist <- mr.deltaG_R(y, X, beta)
    # grad <- -logELCensgrad_R(omegas, deltas, epsilons, lambda, gradlist, weights) # negative gradient
    # attr(res, "gradient") <- grad
    # attr(res, "gradient") <- grad(-logEL.smooth_R,)
  }
  else {
    # TODO: if not converged, what should be the gradient...??
    res <- Inf
    attr(res, "gradient") <- Inf
  }
  return(res)
}