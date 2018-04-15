#' Posterior sampler for mean regression (location model)
#' 
#' @param y Length-\code{nObs} vector of response values.
#' @param X \code{nObs x nEqs} matrix of constraints.
#' @param nsamples Number of samples to obtain.
#' @param nburn number of samples to discard before saving the chain.
#' @param betaInit Length-\code{nEqs} vector of initial value for the chain. 
#' @param sigs Length-\code{nObs} vector of tuning parameters. 
#' @param max_iter Maximum number of Newton-Raphson steps.
#' @param rel_tol Relative tolerance of Newton-Raphson convergence.
#' @return \code{nEqs x nsamples} matrix of Markov Chain.
#' @details ...
#' @export mr.post
mr.post <- function(y, X, nsamples, nburn, betaInit, sigs, 
                    max_iter = 100, rel_tol = 1e-7) {
  # input checks
  if(nrow(X) != length(y)) {
    stop("X and y have inconsistent dimensions.")
  }
  if(ncol(X) != length(betaInit)) {
    stop("X and beta have inconsistent dimensions.")
  }
  .MeanReg_post(y, t(X), nsamples, nburn, betaInit, sigs, 
                 maxIter = max_iter, relTol = rel_tol)
}