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
omega.hat <- function(G, deltas, epsilons, max_iter = 100, 
                      rel_tol = 1e-7, abs_tol = 1e-3, support = FALSE, verbose = FALSE) {
  if (missing(deltas) && missing(epsilons)) {
    lambda <- .lambdaNR(t(G), maxIter = max_iter, relTol = rel_tol, 
                        support = support, verbose = verbose)
    omegahat <- .omega.hat(t(G), lambda = lambda, support = support)
  }
  else {
    # Note: inital omegas obtained from non-censored optimization
    lambda <- .lambdaNR(t(G), maxIter = 100, relTol = rel_tol, 
                        support = support, verbose = verbose)
    omegasInit <- .omega.hat(t(G), lambda = lambda, support = support)
    if (any(is.nan(omegasInit))) {
      # stop("Initial omegas are nans.")
      return(rep(NaN,length(deltas)))
    }
    omegahat <- .omega.hat.EM(omegasInit, t(G), deltas, epsilons,
                              max_iter, rel_tol, absTol = abs_tol,
                              support = support, verbose = verbose)
  }
  return(omegahat)
}

# omega.hat <- function(G, lambda, deltas, epsilons) {
#   if (missing(deltas) && missing(epsilons)) {
#     omegahat <- .omega.hat(t(G), lambda)
#   }
#   else {
#     ...
#   }
# }