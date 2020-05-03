#' Returns the empirical distribution, omegas, from empirical likelihood (EL) maximization using an EM algorithm with discontinuitu correction.
#'
#' @param G \code{n_obs x n_eqs} matrix of constraints.
#' @param deltas Length-\code{n_obs} vector of censor indicators, omit if no censoring.
#' @param epsilons Length-\code{n_obs} vector of residuals, omit if no censoring.
#' @param max_iter Maximum number of Newton-Raphson steps.
#' @param abs_tol Absolute tolerance of Newton-Raphson convergence.
#' @param verbose Display number of steps and tolerance criterion when algorithm terminates.
#' @return Length-\code{n_eqs} vector for the resulting empirical distribution, omegas, if the optimization algorithm converged; 1/n_obs if did not converge.
#' @example examples/omega_hat_EM_smooth.R
#' @export omega_hat_EM_smooth
omega_hat_EM_smooth <- function(G, deltas, epsilons, sp=10, max_iter = 100, 
                                rel_tol = 1e-5, abs_tol = 1e-3, support = FALSE, verbose = FALSE) {
  # Note: inital omegas obtained from non-censored optimization
  lambda <- .LambdaNR(t(G), max_iter = max_iter, rel_tol = rel_tol, support = support, verbose = verbose)
  omegas_init <- .OmegaHat(t(G), lambda, support)
  if (any(is.nan(omegas_init))) {
    # stop("Initial omegas are nans.")
    return(rep(NaN,length(deltas)))
  }
  omegahat <- .OmegaHatEMSmooth(omegas_init, t(G), deltas, epsilons,
                                sp, max_iter, rel_tol, abs_tol, support, verbose)
  return(omegahat)
}