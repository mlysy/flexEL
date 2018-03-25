#' log empirical likelihood for quantile regression
#'
#' @param X \code{nVars x nObs} matrix of constraints.
#' @param y Length-\code{nObs} vector of response values.
#' @param beta Length-\code{nVars} vector of coefficients. 
#' @param alpha Quantile level, a calar between 0 and 1. 
#' @param max_iter Maximum number of Newton-Raphson steps.
#' @param rel_tol Relative tolerance of Newton-Raphson convergence.
#' @return Log empirical likelihood of the input beta
#' @details ...
#' @export
qr.logel <- function(y, X, alpha, beta, max_iter = 100, rel_tol = 1e-7) {
      # input checks
      if(nrow(X) != length(y)) {
            stop("X and y have inconsistent dimensions.")
      }
      if(ncol(X) != length(beta)) {
            stop("X and beta have inconsistent dimensions.")
      }
      if(alpha <= 0 || alpha >= 1) {
            stop("alpha should be between 0 and 1.")
      }
      .QuantReg_logEL(y, t(X), alpha, beta, maxIter = max_iter, relTol = rel_tol)
}
