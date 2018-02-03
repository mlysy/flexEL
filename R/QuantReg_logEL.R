#' log empirical likelihood for quantile regression
#'
#' @param nObs Number of observations.
#' @param nEqs Number of estimating equations.  
#' @param X \code{nVars x nObs} matrix of constraints.
#' @param y Length-\code{nObs} vector of response values.
#' @param beta Length-\code{nVars} vector of coefficients. 
#' @param alpha Quantile level, a calar between 0 and 1. 
#' @param lambda0 Length-\code{nEqs} vector of initial value of lambda.
#' @param max_iter Maximum number of Newton-Raphson steps.
#' @param eps Relative tolerance of Newton-Raphson convergence.
#' @return Log empirical likelihood of the input beta
#' @details ...
#' @export
logEL_QuantReg <- function(X, y, nObs, nEqs, beta, alpha, lambda0, max_iter = 100, eps = 1e-7) {
      # input checks
      if(ncol(X) != length(y)) {
            stop("X and y have inconsistent dimensions.")
      }
      if(nrow(X) != length(beta)) {
            stop("X and beta have inconsistent dimensions.")
      }
      if(alpha <= 0 || alpha >= 1) {
            stop("alpha should be between 0 and 1.")
      }
      .logEL_QuantReg(y, X, nObs, nEqs, beta, alpha, lambda0, maxIter = 100L, eps = 1e-7)
}
