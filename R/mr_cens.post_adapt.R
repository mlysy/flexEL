#' Posterior sampler for mean regression (location model)
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
#' @export mr_cens.post_adapt
mr_cens.post_adapt <- function(y, X, deltas, nsamples, nburn, betaInit, 
                               mwgSd, rvDoMcmc, max_iter = 100, rel_tol = 1e-7) {
  # # input conversion
  # if (is.vector(BetaInit)) BetaInit <- matrix(BetaInit,length(BetaInit),1)
  # if (is.vector(mwgSd)) mwgSd <- matrix(mwgSd,length(mwgSd),1)
  # input checks
  if(nrow(X) != length(y)) {
    stop("X and y have inconsistent dimensions.")
  }
  if(ncol(X) != length(betaInit)) {
    stop("X and beta have inconsistent dimensions.")
  }
  # nBet <- nrow(BetaInit)
  # numBet <- ncol(BetaInit)
  
  # obtain initial value for first EM from uncensored case
  G <- .MeanReg_evalG(y,t(X), betaInit)
  lambda <- .lambdaNR(G, maxIter = max_iter, relTol = rel_tol, verbose = FALSE)
  omegasInit <- .omega.hat(G, lambda)
  if (anyNA(omegasInit)) {
    stop("Initial omegas are nans.")
  }
  if (missing(rvDoMcmc)) rvDoMcmc <- rep(1,length(betaInit))
  .MeanRegCens_post_adapt(omegasInit, y, t(X), deltas, nsamples, nburn, betaInit, 
                          mwgSd, rvDoMcmc, maxIter = max_iter, relTol = rel_tol)
}