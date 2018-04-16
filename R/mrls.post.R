#' Posterior sampler for mean regression (location-scale model)
#' 
#' @param y Length-\code{nObs} vector of response values.
#' @param X \code{nObs x nBet} matrix of location function covariates.
#' @param Z \code{nObs x nGam} matrix of scale function covariates.
#' @param nsamples Number of samples to obtain.
#' @param nburn number of samples to discard before saving the chain.
#' @param betaInit Length-\code{nBet} vector of initial value for the chain. 
#' @param gammaInit Length-\code{nGam} vector of initial value for the chain. 
#' @param sigs Length-\code{nEqs} vector of tuning parameters. 
#' @param betalen A scalar of the length of the location model parameter.
#' @param max_iter Maximum number of Newton-Raphson steps.
#' @param rel_tol Relative tolerance of Newton-Raphson convergence.
#' @return \code{nEqs x nsamples} matrix of Markov Chain.
#' @details ...
#' @export mrls.post
mrls.post <- function(y, X, Z, nsamples, nburn, betaInit, gammaInit, 
                      sigs, max_iter = 100, rel_tol = 1e-7) {
  # input checks
  if(nrow(X) != length(y)) {
    stop("X and y have inconsistent dimensions.")
  }
  if(ncol(X) != length(betaInit)) {
    stop("X and beta have inconsistent dimensions.")
  }
  if(nrow(Z) != length(y)) {
    stop("Z and y have inconsistent dimensions.")
  }
  if(ncol(Z) != length(gammaInit)) {
    stop("Z and gamma have inconsistent dimensions.")
  }
  .MeanRegLS_post(y, t(X), t(Z), nsamples, nburn, 
                  betaInit, gammaInit, sigs, 
                  maxIter = max_iter, relTol = rel_tol)
}