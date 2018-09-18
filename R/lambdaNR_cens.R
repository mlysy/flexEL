#' Solve the inner loop optimization of a censored EL function
#'
#' @template args-G
#' @param weights Length-\code{nObs} vector of weights for weighted Newton-Raphson algorithm.
#' @template args-max_iter
#' @template args-rel_tol
#' @template args-verbose
#' @return Length-\code{nEq} vector corresponding to the solution of the optimization problem.
#' @details The inner-loop optimization of EL is ...
#' @export lambdaNR_cens
lambdaNR_cens <- function(G, weights, max_iter = 100, rel_tol = 1e-7, verbose = FALSE) { 
  lambda <- .lambdaNRC(G = t(G), weights, maxIter = max_iter, relTol = rel_tol, verbose = verbose)
  return(lambda)
}

