#' Solves the inner optimization problem of an EL maximization.
#'
#' @template args-G
#' @template args-max_iter
#' @template args-rel_tol
#' @template args-verbose
#' @example examples/lambdaNR.R
#' @return Length-\code{nEq} vector corresponding to the solution of the optimization problem.
#' @export lambdaNR
lambdaNR <- function(G, max_iter = 100, rel_tol = 1e-7, support = FALSE, verbose = FALSE) { 
  lambda <- .LambdaNR(G = t(G), 
                      max_iter = max_iter, rel_tol = rel_tol, support = support, 
                      verbose = verbose)
  return(lambda)
}

# LambdaNR <- function(G, weights, max_iter = 100, rel_tol = 1e-7, verbose = FALSE) { 
#     # check whether weights is given and call the corresponding NR funciton
#     if (missing(weights)) {
#         ans <- .LambdaNR(G = t(G), max_iter = max_iter, rel_tol = rel_tol, verbose = verbose)
#     }
#     else {
#         ans <- .LambdaNRC(G = t(G), weights, max_iter = max_iter, rel_tol = rel_tol, verbose = verbose)
#     }
#     # check convergence of NR
#     if(ans$convergence) {
#         ans <- ans$lambda
#     } 
#     else {
#         ans <- rep(NaN, ncol(G)) # ncol(G) == nEqs
#     }
#     ans
# }
