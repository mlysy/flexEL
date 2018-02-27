#' Returns the omegas, i.e., the empirical probabilities.
#'
#' @param G \code{nObs x nEqns} matrix of constraints.
#' @param max_iter Maximum number of Newton-Raphson steps.
#' @param eps Relative tolerance of Newton-Raphson convergence.
#' @param verbose Display number of steps and tolerance criterion when algorithm terminates.
#' @return Length-\code{nEq} vector corresponding to the solution of the optimization problem.
#' @details The inner-loop optimization of EL is ...
#' @export
omega.hat <- function(G, deltas, epsilons, max_iter = 100, rel_tol = 1e-7, verbose = FALSE) {
    if (missing(deltas) && missing(epsilons)) {
        omega.hat <- .omega.hat(t(G), max_iter, rel_tol, verbose)
    }
    else {
        omega.hat <- .omega.hat.EM(...)
    }
    omega.hat
}