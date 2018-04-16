#' Posterior sampler for quantile regression (location model)
#' 
#' @param y Length-\code{nObs} vector of response values.
#' @param X \code{nObs x nEqs} matrix of constraints.
#' @param alphas a vector of quantile levels.
#' @param nsamples Number of samples to obtain.
#' @param nburn number of samples to discard before saving the chain.
#' @param betaInit Length-\code{nEqs} vector of initial value for the chain. 
#' @param sigs Length-\code{nEqs} vector of tuning parameters. 
#' @param max_iter Maximum number of Newton-Raphson steps.
#' @param rel_tol Relative tolerance of Newton-Raphson convergence.
#' @return \code{nEqs x nsamples} matrix of Markov Chain.
#' @details ...
#' @export qr.post
qr.post <- function(y, X, alphas, nsamples, nburn, betaInit, sigs, 
                    max_iter = 100, rel_tol = 1e-7) {
  # input checks
  if(nrow(X) != length(y)) {
      stop("X and y have inconsistent dimensions.")
  }
  if(ncol(X) != length(betaInit)) {
      stop("X and beta have inconsistent dimensions.")
  }
  # add a first entry of alpha as the number of quantile levels
  alpha <- c(length(alphas), alphas) 
  .QuantReg_post(y, t(X), alpha, nsamples, nburn, betaInit, sigs, 
                 maxIter = max_iter, relTol = rel_tol)
}