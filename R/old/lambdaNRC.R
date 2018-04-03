#' Solve the inner loop optimization of an EL function (censoring).
#'
#' @param G \code{nObs x nEqns} matrix of constraints.
#' @param max_iter Maximum number of Newton-Raphson steps.
#' @param eps Relative tolerance of Newton-Raphson convergence.
#' @param verbose Display number of steps and tolerance criterion when algorithm terminates.
#' @return Length-\code{nEq} vector corresponding to the solution of the optimization problem.
#' @details The inner-loop optimization of EL is ...
#' @export
lambdaNRC <- function(G, qs, max_iter = 100, eps = 1e-7, verbose = FALSE) {
    nObs <- nrow(G)
    nEqs <- ncol(G)
    ans <- .lambdaNRC(G = t(G), qs, maxIter = max_iter, eps = eps, verbose = verbose)
    if(ans$convergence) {
        ans <- ans$lambda
    } else {
        ans <- rep(NA, nEqs)
    }
    ans
}
