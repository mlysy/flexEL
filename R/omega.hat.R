#' Returns the empirical distribution, omegas, from empirical likelihood (EL) maximization.
#'
#' @param G \code{nObs x nEqns} matrix of constraints.
#' @param deltas Length-\code{nObs} vector of censor indicators, omit if no censoring.
#' @param epsilons Length-\code{nObs} vector of residuals, omit if no censoring.
#' @param max_iter Maximum number of Newton-Raphson steps.
#' @param rel_tol Relative tolerance of Newton-Raphson convergence.
#' @param verbose Display number of steps and tolerance criterion when algorithm terminates.
#' @return Length-\code{nEq} vector for the resulting empirical distribution, omegas, if the optimization algorithm converged; 1/nObs if did not converge.
#' @details The inner-loop optimization of EL is ...
#' @export
omega.hat <- function(G, deltas, epsilons, max_iter = 100, rel_tol = 1e-7, verbose = FALSE) {
    if (missing(deltas) && missing(epsilons)) {
        omegaOut <- .omega.hat(t(G), max_iter, rel_tol, verbose)
    }
    else {
        omegaOut <- .omega.hat.EM(t(G), deltas, epsilons, 
                                  max_iter, rel_tol, verbose)
    }
    if (omegaOut$convergence) {
        return(omegaOut$omegas)
    }
    else { 
        # TODO: should give warning here, verbose???
        n <- nrow(G)
        return(rep(1/n,n))
    }
}