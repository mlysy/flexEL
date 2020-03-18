#' Returns the empirical distribution, omegas, from empirical likelihood (EL) maximization.
#'
#' @param G \code{n_obs x n_eqs} matrix of constraints.
#' @param deltas Length-\code{n_obs} vector of censor indicators, omit if no censoring.
#' @param epsilons Length-\code{n_obs} vector of residuals, omit if no censoring.
#' @param max_iter Maximum number of Newton-Raphson steps.
#' @param rel_tol Relative tolerance of Newton-Raphson convergence.
#' @param verbose Display number of steps and tolerance criterion when algorithm terminates.
#' @return Length-\code{n_eqs} vector for the resulting empirical distribution, omegas, if the optimization algorithm converged; 1/n_obs if did not converge.
#' @details The inner-loop optimization of EL is ...
#' @export omega_hat
omega_hat <- function(G, deltas, epsilons, 
                      max_iter = 100, rel_tol = 1e-7, abs_tol = 1e-3, support = FALSE, verbose = FALSE) {
  if (missing(deltas) && missing(epsilons)) {
    lambda <- .LambdaNR(t(G), 
                        max_iter = max_iter, 
                        rel_tol = rel_tol, 
                        support = support, 
                        verbose = verbose)
    omegahat <- .OmegaHat(t(G), lambda = lambda, support = support)
  }
  else {
    # Note: inital omegas obtained from non-censored optimization
    lambda <- .LambdaNR(t(G), 
                        max_iter = max_iter, 
                        rel_tol = rel_tol, 
                        support = support, 
                        verbose = verbose)
    omegas_init <- .OmegaHat(t(G), lambda = lambda, support = support)
    if (any(is.nan(omegas_init))) {
      # stop("Initial omegas are nans.")
      return(rep(NaN,length(deltas)))
    }
    omegahat <- .OmegaHatEM(omegas_init, t(G), deltas, epsilons,
                            max_iter, rel_tol, abs_tol = abs_tol,
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