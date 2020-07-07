#' Solves the inner optimization problem of an EL maximization.
#'
#' @template args-G
#' @template args-max_iter
#' @template args-rel_tol
#' @template args-verbose
#' @example examples/lambdaNR.R
#' @return Vector of length `n_eq` corresponding to the solution of the optimization problem.
#' @export lambdaNR
lambdaNR <- function(G, max_iter = 100, rel_tol = 1e-7, support = FALSE, verbose = FALSE) {
  lambda <- .LambdaNR(G = t(G),
                      max_iter = max_iter, rel_tol = rel_tol, support = support,
                      verbose = verbose)
  return(lambda)
}
