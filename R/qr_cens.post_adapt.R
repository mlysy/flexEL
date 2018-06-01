#' Posterior sampler for quantile regression with censoring (location model)
#' 
#' @param y Length-\code{nObs} vector of response values.
#' @param X \code{nObs x nEqs} matrix of constraints.
#' @param deltas Length-\code{nObs} vector of censoring indicators.
#' @param alpha a scalar of quantile level (TODO: multiple quantile levels).
#' @param nsamples Number of samples to obtain.
#' @param nburn number of samples to discard before saving the chain.
#' @param BetaInit \code{nEqs x numBeta} matrix of initial value for the chain. Each column is a beta vector. If a vector is passed to this function, it will be converted to a matrix with only one column.
#' @param Sigs \code{nEqs x numBeta} matrix of tuning parameters. 
#' @param rvDoMcmc \code{nEqs x numBeta} matrix of indicators (1 or 0) for whether or not updating the corresponding parameter.
#' @param max_iter Maximum number of Newton-Raphson steps.
#' @param rel_tol Relative tolerance of Newton-Raphson convergence.
#' @return \code{nEqs x nsamples} matrix of Markov Chain.
#' @details ...
#' @export qr_cens.post_adapt
qr_cens.post_adapt <- function(y, X, deltas, alpha, nsamples, nburn, BetaInit, 
                               Sigs, RvDoMcmc, DoAdapt,
                               max_iter = 100, rel_tol = 1e-7) {
  # # input conversion
  if (is.vector(BetaInit)) BetaInit <- matrix(BetaInit,length(BetaInit),1)
  if (is.vector(Sigs)) Sigs <- matrix(Sigs,length(Sigs),1)
  if (missing(RvDoMcmc)) {
    RvDoMcmc <- matrix(1, nrow = nrow(BetaInit), ncol = ncol(BetaInit))
  }
  if (missing(DoAdapt)) {
    DoAdapt <- matrix(1, nrow = nrow(BetaInit), ncol = ncol(BetaInit))
  }
  if (is.vector(RvDoMcmc)) RvDoMcmc <- matrix(RvDoMcmc, length(RvDoMcmc), 1)
  if (is.vector(DoAdapt)) DoAdapt <- matrix(RvDoMcmc, length(RvDoMcmc), 1)
  # input checks
  if(nrow(X) != length(y)) {
    stop("X and y have inconsistent dimensions.")
  }
  if(ncol(X) != length(betaInit)) {
    stop("X and beta have inconsistent dimensions.")
  }
  # if(!all(dim(BetaInit) == dim(Sigs))) {
  #   stop("BetaInit and Sigs have inconsistent dimensions.")
  # }
  # if(ncol(X) != nrow(BetaInit)) {
  #   stop("X and beta have inconsistent dimensions.")
  # }
  
  # obtain initial value for first EM from uncensored case
  G <- .QuantReg_evalG(y, t(X), c(1,alpha), BetaInit)
  lambda <- .lambdaNR(G, maxIter = max_iter, relTol = rel_tol, verbose = FALSE)
  omegasInit <- .omega.hat(G, lambda)
  # if (missing(rvDoMcmc)) rvDoMcmc <- rep(1,length(BetaInit))
  .QuantRegCens_post_adapt(omegasInit, y, t(X), deltas, c(1,alpha), nsamples, nburn, 
                           BetaInit, Sigs, RvDoMcmc, DoAdapt, 
                           maxIter = max_iter, relTol = rel_tol)
}