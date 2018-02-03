#' log empirical likelihood for mean regression with location scale model
#'
#' @param nObs Number of observations.
#' @param nEqs Number of estimating equations.
#' @param X \code{nVars x nObs} matrix of constraints.
#' @param y Length-\code{nObs} vector of response values.
#' @param theta Length-\code{nVars} vector of coefficients.
#' @param lambda0 Length-\code{nEqs} vector of initial value of lambda.
#' @param max_iter Maximum number of Newton-Raphson steps.
#' @param eps Relative tolerance of Newton-Raphson convergence.
#' @return Log empirical likelihood of the input theta
#' @details ...
#' @export
mrls.logel <- function(y, X, theta, max_iter = 100, eps = 1e-7) {
    # input checks
    if(ncol(X) != length(y)) {
        stop("X and y have inconsistent dimensions.")
    }
    if(nrow(X) != length(theta)) {
        stop("X and theta have inconsistent dimensions.")
    }
    .MeanRegLS_logEL(y, X, theta, maxIter = 100L, eps = 1e-7)
}
