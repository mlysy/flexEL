#' Solve the inner loop optimization of an EL function
#'
#' @template args-G
#' @template args-max_iter
#' @template args-rel_tol
#' @template args-verbose
#' @return Length-\code{nEq} vector corresponding to the solution of the optimization problem.
#' @details The inner-loop optimization of EL is ...
#' @export lambdaNR
lambdaNR <- function(G, max_iter = 100, rel_tol = 1e-7, verbose = FALSE) { 
  lambda <- .lambdaNR(G = t(G), maxIter = max_iter, relTol = rel_tol, verbose = verbose)
  return(lambda)
}

# lambdaNR <- function(G, weights, max_iter = 100, rel_tol = 1e-7, verbose = FALSE) { 
#     # check whether weights is given and call the corresponding NR funciton
#     if (missing(weights)) {
#         ans <- .lambdaNR(G = t(G), maxIter = max_iter, relTol = rel_tol, verbose = verbose)
#     }
#     else {
#         ans <- .lambdaNRC(G = t(G), weights, maxIter = max_iter, relTol = rel_tol, verbose = verbose)
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
