#' Posterior sampler for mean regression (location-scale model)
#' 
#' @param y Length-\code{nObs} vector of response values.
#' @param X \code{nObs x nBet} matrix of location function covariates.
#' @param Z \code{nObs x nGam} matrix of scale function covariates.
#' @param nsamples Number of samples to obtain.
#' @param nburn number of samples to discard before saving the chain.
#' @param BetaInit \code{nBet x numBeta} matrix of initial value for the chain. 
#' @param GammaInit \code{nGam x numGamma} matrix of initial value for the chain. 
#' @param Sigs Length-\code{nEqs} vector of tuning parameters. 
#' @param max_iter Maximum number of Newton-Raphson steps.
#' @param rel_tol Relative tolerance of Newton-Raphson convergence.
#' @return \code{nEqs x nsamples} matrix of Markov Chain.
#' @details ...
#' @export mrls.post
mrls.post <- function(y, X, Z, nsamples, nburn, BetaInit, GammaInit, 
                      Sigs, max_iter = 100, rel_tol = 1e-7) {
  # input conversion
  if (is.vector(BetaInit)) BetaInit <- matrix(BetaInit,length(BetaInit),1)
  if (is.vector(GammaInit)) GammaInit <- matrix(GammaInit,length(GammaInit),1)
  if (is.vector(Sigs)) Sigs <- matrix(Sigs,length(Sigs),1)
  # input checks
  if(nrow(X) != length(y)) {
    stop("X and y have inconsistent dimensions.")
  }
  if(ncol(X) != nrow(BetaInit)) {
    stop("X and beta have inconsistent dimensions.")
  }
  if(nrow(Z) != length(y)) {
    stop("Z and y have inconsistent dimensions.")
  }
  if(ncol(Z) != nrow(GammaInit)) {
    stop("Z and gamma have inconsistent dimensions.")
  }
  .MeanRegLS_post(y, t(X), t(Z), nsamples, nburn, 
                  BetaInit, GammaInit, Sigs, 
                  maxIter = max_iter, relTol = rel_tol)
}