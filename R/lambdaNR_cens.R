#' Solve the inner loop optimization of a censored EL function
#'
#' @template args-G
#' @param weights Length-\code{n_obs} vector of weights for weighted Newton-Raphson algorithm.
#' @template args-max_iter
#' @template args-rel_tol
#' @template args-verbose
#' @return Length-\code{n_eqs} vector corresponding to the solution of the optimization problem.
#' @details The inner-loop optimization of EL is ...
#' @export lambdaNR_cens
lambdaNR_cens <- function(G, weights, max_iter = 100, rel_tol = 1e-7, support = FALSE,
                          verbose = FALSE) { 
  lambda <- .LambdaNRC(G = t(G), weights, max_iter = max_iter, rel_tol = rel_tol, 
                       support = support, verbose = verbose)
  return(lambda)
}

