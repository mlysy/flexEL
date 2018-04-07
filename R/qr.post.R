#' Posterior sampler for quantile regression (location model)
#' 
#' @param y Length-\code{nObs} vector of response values.
#' @param X \code{nObs x nEqs} matrix of constraints.
#' @param beta Length-\code{nVars} vector of coefficients.
#' @param max_iter Maximum number of Newton-Raphson steps.
#' @param rel_tol Relative tolerance of Newton-Raphson convergence.
#' @return Log empirical likelihood of the input beta
#' @details ...
#' @export
qr.post <- function(y, X, alpha, nsamples, nburn, betaInit, sigs, 
                    max_iter = 100, rel_tol = 1e-7) {
    # input checks
    if(nrow(X) != length(y)) {
        stop("X and y have inconsistent dimensions.")
    }
    if(ncol(X) != length(betaInit)) {
        stop("X and beta have inconsistent dimensions.")
    }
    .QuantReg_post(y, t(X), alpha, nsamples, nburn, betaInit, sigs, 
                   maxIter = max_iter, relTol = rel_tol)
}