#' Solve the inner loop optimization of an EL function.
#'
#' @param G \code{nObs x nEqns} matrix of constraints.
#' @param weights Optional length-\code{nObs} vector of weights for weighted Newton-Raphson algorithm.
#' @param max_iter Maximum number of Newton-Raphson steps.
#' @param rel_tol Relative tolerance of Newton-Raphson convergence.
#' @param verbose Display number of steps and tolerance criterion when algorithm terminates.
#' @return Length-\code{nEq} vector corresponding to the solution of the optimization problem.
#' @details The inner-loop optimization of EL is ...
#' @export
lambdaNR <- function(G, weights, max_iter = 100, rel_tol = 1e-7, verbose = FALSE) { 
    # check whether weights is given and call the corresponding NR funciton
    if (missing(weights)) {
        ans <- .lambdaNR(G = t(G), maxIter = max_iter, relTol = rel_tol, verbose = verbose)
    }
    else {
        ans <- .lambdaNRC(G = t(G), weights, maxIter = max_iter, relTol = rel_tol, verbose = verbose)
    }
    # check convergence of NR
    if(ans$convergence) {
        ans <- ans$lambda
    } 
    else {
        ans <- rep(NaN, ncol(G)) # ncol(G) == nEqs
    }
    ans
}