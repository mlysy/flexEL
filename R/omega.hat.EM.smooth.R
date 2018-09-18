#' Returns the empirical distribution, omegas, from empirical likelihood (EL) maximization.
#'
#' @param G \code{nObs x nEqns} matrix of constraints.
#' @param deltas Length-\code{nObs} vector of censor indicators, omit if no censoring.
#' @param epsilons Length-\code{nObs} vector of residuals, omit if no censoring.
#' @param max_iter Maximum number of Newton-Raphson steps.
#' @param abs_tol Absolute tolerance of Newton-Raphson convergence.
#' @param verbose Display number of steps and tolerance criterion when algorithm terminates.
#' @return Length-\code{nEq} vector for the resulting empirical distribution, omegas, if the optimization algorithm converged; 1/nObs if did not converge.
#' @details The inner-loop optimization of EL is ...
#' @export omega.hat.EM.smooth
omega.hat.EM.smooth <- function(G, deltas, epsilons, sp=10, max_iter = 100, 
                                rel_tol = 1e-5, abs_tol = 1e-3, verbose = FALSE) {
  # Note: inital omegas obtained from non-censored optimization
  lambda <- .lambdaNR(t(G), maxIter = 100, relTol = rel_tol, verbose = verbose)
  omegasInit <- .omega.hat(t(G), lambda)
  if (any(is.nan(omegasInit))) {
    # stop("Initial omegas are nans.")
    return(rep(NaN,length(deltas)))
  }
  omegahat <- .omega.hat.EM.smooth(omegasInit, t(G), deltas, epsilons,
                                   sp, max_iter, rel_tol, abs_tol, verbose)
  return(omegahat)
}