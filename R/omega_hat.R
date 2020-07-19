#' Returns the empirical distribution, omegas, from empirical likelihood (EL) maximization.
#'
#' @template arg-G
#' @template arg-delta
#' @template arg-eps
#' @template args-lambda_precision
#' @template arg-abs_tol
#' @param verbose Display number of steps and tolerance criterion when algorithm terminates.
#' @example examples/omega_hat.R
#' @return Length-`n_eqs` vector for the resulting empirical distribution, omegas, if the optimization algorithm converged; 1/n_obs if did not converge.
#' @noRd
omega_hat <- function(G, delta = NULL, eps = NULL, support = FALSE, 
                      max_iter = 100, rel_tol = 1e-7, abs_tol = 1e-3, 
                      verbose = FALSE) {
  
  if (is.null(delta) && is.null(eps)) {
    lambda <- .LambdaNR(t(G), 
                        max_iter = max_iter, 
                        rel_tol = rel_tol, 
                        support = support, 
                        verbose = verbose)
    return(.OmegaHat(t(G), lambda = lambda, support = support))
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
      return(rep(NaN, length(delta)))
    }
    else {
      return(.OmegaHatEM(omegas_init, t(G), delta, eps,
                         max_iter, rel_tol, abs_tol = abs_tol,
                         support = support, verbose = verbose))
    }
  }
}
