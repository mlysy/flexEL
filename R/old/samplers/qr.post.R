#' Posterior sampler for quantile regression (location model)
#' 
#' @param y Length-\code{nObs} vector of response values.
#' @param X \code{nObs x nEqs} matrix of constraints.
#' @param alphas a vector of quantile levels.
#' @param nsamples Number of samples to obtain.
#' @param nburn number of samples to discard before saving the chain.
#' @param BetaInit \code{nEqs x numBeta} matrix of initial value for the chain. Each column is a beta vector. If a vector is passed to this function, it will be converted to a matrix with only one column.
#' @param Sigs \code{nEqs x numBeta} matrix of tuning parameters. 
#' @param RvDoMcmc \code{nEqs x numBeta} matrix of indicators (1 or 0) for whether or not updating the corresponding parameter.
#' @param max_iter Maximum number of Newton-Raphson steps.
#' @param rel_tol Relative tolerance of Newton-Raphson convergence.
#' @return \code{nEqs x nsamples} matrix of Markov Chain.
#' @details ...
#' @export qr.post
qr.post <- function(y, X, alphas, nsamples, nburn, BetaInit, Sigs, RvDoMcmc, 
                    max_iter = 100, rel_tol = 1e-7) {
  # input conversion
  if (is.vector(BetaInit)) BetaInit <- matrix(BetaInit,length(BetaInit),1)
  if (is.vector(Sigs)) Sigs <- matrix(Sigs,length(Sigs),1)
  # if RvDoMcmc is not specified, then all get updated
  if (missing(RvDoMcmc)) {
    RvDoMcmc <- matrix(1, nrow = nrow(BetaInit), ncol = ncol(BetaInit))
  }
  if (is.vector(RvDoMcmc)) RvDoMcmc <- matrix(RvDoMcmc, length(RvDoMcmc), 1)
  # add a first entry of alpha as the number of quantile levels
  alpha <- c(length(alphas), alphas) 
  # input checks
  if(nrow(X) != length(y)) {
    stop("X and y have inconsistent dimensions.")
  }
  if(!all(dim(BetaInit) == dim(Sigs))) {
    stop("BetaInit and Sigs have inconsistent dimensions.")
  }
  if(!all(dim(BetaInit) == dim(RvDoMcmc))) {
    stop("BetaInit and RvDoMcmc have inconsistent dimensions.")
  }
  if(ncol(X) != nrow(BetaInit)) {
    stop("X and beta have inconsistent dimensions.")
  }
  .QuantReg_post(y, t(X), alpha, nsamples, nburn, BetaInit, Sigs, RvDoMcmc, 
                 maxIter = max_iter, relTol = rel_tol)
}