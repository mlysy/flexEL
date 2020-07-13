#' Solves the inner optimization problem of an EL maximization problem.
#'
#' @template arg-G
#' @param weights Optional length-`n_obs` vector of weights for weighted Newton-Raphson algorithm for right-censored data.
#' @template args-lambda_precision
#' @template arg-support
#' @template arg-verbose
#' @example examples/lambdaNR.R
#' @return Vector of length `n_eq` corresponding to the solution of the optimization problem.
#' @export lambdaNR
lambdaNR <- function(G, weights = NULL, support = FALSE, max_iter = 100, rel_tol = 1e-7, verbose = FALSE) {
  if (is.null(weights)) {
    lambda <- .LambdaNR(G = t(G),
                        max_iter = max_iter, rel_tol = rel_tol, 
                        support = support, verbose = verbose)
  }
  else {
    lambda <- .LambdaNRCens(G = t(G), weights = weights, 
                            max_iter = max_iter, rel_tol = rel_tol, 
                            support = support, verbose = verbose)
  }
  return(lambda)
}
