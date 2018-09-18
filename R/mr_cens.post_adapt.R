#' Posterior sampler for mean regression (location model)
#' 
#' @template args-y
#' @template args-X
#' @template args-delta
#' @template args-nsamp
#' @template args-nburn
#' @template args-betaInit
#' @template args-mwgSd
#' @template args-doMcmc
#' @template args-max_iter
#' @template args-rel_tol
#' @return \code{nEqs x nsamp} matrix of posterior samples.
#' @details ...
#' @export mr_cens.post_adapt
mr_cens.post_adapt <- function(y, X, delta, nsamp, nburn, betaInit, 
                               mwgSd, doMcmc, DoAdapt, 
                               max_iter = 100, rel_tol = 1e-7, abs_tol = 1e-3) {
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
  if (missing(doMcmc)) doMcmc <- rep(1,length(betaInit))
  if (missing(DoAdapt)) DoAdapt <- rep(1,length(betaInit))
  .MeanRegCens_post_adapt(omegasInit, y, t(X), delta, nsamp, nburn, betaInit, 
                          mwgSd, doMcmc, DoAdapt, 
                          maxIter = max_iter, relTol = rel_tol, absTol = abs_tol)
}