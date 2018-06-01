#' Adaptive posterior sampler for qunatile regression (location-scale model)
#' 
#' @param y Length-\code{nObs} vector of response values.
#' @param X \code{nObs x nBet} matrix of location function covariates.
#' @param Z \code{nObs x nGam} matrix of scale function covariates.
#' @param alpha a scalar of quantile level (TODO: multiple quantile levels).
#' @param nsamples Number of samples to obtain.
#' @param nburn number of samples to discard before saving the chain.
#' @param betaInit \code{nBet} vector of initial value for the parameters in the linear location function. 
#' @param gammaInit \code{nGam} vector of initial value for the parameters in the exponential scale function. 
#' @param sig2Init A scalar of initial value for scale parameter.
#' @param nuInit A scalar of initial value for the quantile value of the error term.
#' @param mwgSd Length-\code{nEqs} vector of tuning parameters. 
#' @param rvDoMcmc \code{nEqs x numBeta} matrix of indicators (1 or 0) for whether or not updating the corresponding parameter.
#' @param max_iter Maximum number of Newton-Raphson steps.
#' @param rel_tol Relative tolerance of Newton-Raphson convergence.
#' @return \code{nEqs x nsamples} matrix of Markov Chain.
#' @details ...
#' @export qrls_cens.post_adapt
qrls_cens.post_adapt <- function(y, X, Z, deltas, alpha, nsamples, nburn, 
                                 betaInit, gammaInit, sig2Init, nuInit, 
                                 mwgSd, rvDoMcmc, doAdapt,
                                 max_iter = 100, rel_tol = 1e-7) {
  # input conversion
  # if (is.vector(BetaInit)) BetaInit <- matrix(BetaInit,length(BetaInit),1)
  # if (is.vector(GammaInit)) GammaInit <- matrix(GammaInit,length(GammaInit),1)
  # if (is.vector(mwgSd)) mwgSd <- matrix(mwgSd,length(mwgSd),1)
  # if RvDoMcmc is not specified, then all get updated
  # if (missing(RvDoMcmc)) {
  #   RvDoMcmc <- matrix(1, nrow = nrow(BetaInit)+nrow(GammaInit)+1,
  #                      ncol = ncol(BetaInit))
  # }
  # if (is.vector(RvDoMcmc)) RvDoMcmc <- matrix(RvDoMcmc,length(RvDoMcmc),1)
  if (missing(rvDoMcmc)) {
    rvDoMcmc <- rep(1, length(betaInit)+length(gammaInit)+2)
    # RvDoMcmc <- rep(1, nrow(betaInit)+nrow(gammaInit)+2)
  }
  if (missing(doAdapt)) {
    doAdapt <- rep(1, length(betaInit)+length(gammaInit)+2)
    # RvDoMcmc <- rep(1, nrow(betaInit)+nrow(gammaInit)+2)
  }
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
  # obtain initial value for first EM from uncensored case
  G <- .QuantRegLS_evalG(y,t(X),t(Z),betaInit,gammaInit,sig2Init)
  lambda <- .lambdaNR(G, maxIter = max_iter, relTol = rel_tol, verbose = FALSE)
  omegasInit <- .omega.hat(G,lambda)
  if (anyNA(omegasInit)) {
    stop("Initial omegas not converged.")
    # print(omegasInit)
  }
  .QuantRegCensLS_post_adapt(omegasInit, y, t(X), t(Z), deltas, c(1,alpha),
                             nsamples, nburn, 
                             betaInit, gammaInit,sig2Init, nuInit,
                             mwgSd, rvDoMcmc, doAdapt,
                             maxIter = max_iter, relTol = rel_tol)
}