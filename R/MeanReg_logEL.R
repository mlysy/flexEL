#' log empirical likelihood for mean regression
#'
#' @param X \code{nEqs x nObs} matrix of constraints.
#' @param y Length-\code{nObs} vector of response values.
#' @param beta Length-\code{nEqs} vector of coefficients.
#' @param max_iter Maximum number of Newton-Raphson steps.
#' @param rel_tol Relative tolerance of Newton-Raphson convergence.
#' @return Log empirical likelihood of the input beta
#' @details ...
#' @export
mr.logel <- function(y, X, beta, max_iter = 100, rel_tol = 1e-7) {
  # input checks
  if(ncol(X) != length(y)) {
    stop("X and y have inconsistent dimensions.")
  }
  if(nrow(X) != length(beta)) {
    stop("X and beta have inconsistent dimensions.")
  }
  .MeanReg_logEL(y, X, beta, maxIter = max_iter, relTol = rel_tol)
}
