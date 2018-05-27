#' Posterior sampler for mean regression (location-scale model)
#' 
#' @param y Length-\code{nObs} vector of response values.
#' @param X \code{nObs x nEqs} matrix of constraints.
#' @param deltas Length-\code{nObs} vector of censoring indicators.
#' @param nsamples Number of samples to obtain.
#' @param nburn number of samples to discard before saving the chain.
#' @param betaInit \code{nEqs} vector of initial value for the chain.
#' @param mwgSd Length-\code{nObs} vector of tuning parameters. 
#' @param rvDoMcmc Length-\code{nEqs} vector of 1 and 0s, where 1 indicates to update the corresponding entry of beta.
#' @param max_iter Maximum number of Newton-Raphson steps.
#' @param rel_tol Relative tolerance of Newton-Raphson convergence.
#' @return \code{nEqs x nsamples} matrix of Markov Chain.
#' @details ...
#' @export mrls_cens.post_adapt
mrls_cens.post_adapt <- function(y, X, Z, deltas, nsamples, nburn, 
                                 betaInit, gammaInit, sig2Init,
                                 mwgSd, rvDoMcmc, max_iter = 100, rel_tol = 1e-7) {
  # input checks
  if(nrow(X) != length(y)) {
    stop("X and y have inconsistent dimensions.")
  }
  if(ncol(X) != length(betaInit)) {
    stop("X and beta have inconsistent dimensions.")
  }
  # obtain initial value for first EM from uncensored case
  G <- mrls.evalG(y,X,Z,betaInit,gammaInit,sig2Init)
  weights <- evalWeights(deltas,omegas,epsilons)
  lambda <- .lambdaNRC(t(G), weights, maxIter = max_iter, relTol = rel_tol, verbose = FALSE)
  omegasInit <- .omega.hat(t(G), lambda)
  if (anyNA(omegasInit)) {
    stop("Initial omegas not converged.")
    # print(omegasInit)
  }
  if (missing(rvDoMcmc)) rvDoMcmc <- rep(1,length(betaInit))
  .MeanRegCensLS_post_adapt(omegasInit, y, t(X), t(Z), deltas, nsamples, nburn, 
                            betaInit, gammaInit, sig2Init, mwgSd, rvDoMcmc, 
                            maxIter = max_iter, relTol = rel_tol)
}
