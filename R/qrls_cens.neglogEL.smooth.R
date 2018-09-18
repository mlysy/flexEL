#' Calculate the negative log EL with smoothing for quantile regression location-scale model under censoring
#'
#' @param omegas Vector of probability weights
#' @param epsilons Vector of residuals.
#' @param ii An integer indicator between 1 and the length of omegas/epsilons.
#' @return Partial sum of omegas according to residuals that are no larger than the ii-th residual.
#' @details ...
#' @export qrls_cens.neglogEL.smooth
qrls_cens.neglogEL.smooth <- function(y, X, Z, deltas, tau, theta, sp=10) {
  nBet <- ncol(X)
  nGam <- ncol(Z)
  beta <- theta[1:nBet]
  gamma <- theta[(nBet+1):(nBet+nGam)]
  sig2 <- theta[nBet+nGam+1]
  if (sig2 < 0) return(Inf)
  nu <- theta[nBet+nGam+2]
  G <- qrls.evalG.smooth(y, X, Z, tau, beta, gamma, sig2, nu, sp)
  if (anyNA(G)) return(Inf)
  epsilons <- evalEpsilonsLS(y,X,Z,beta,gamma,sig2)
  omegas <- omega.hat.EM.smooth(G,deltas,epsilons,sp)
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