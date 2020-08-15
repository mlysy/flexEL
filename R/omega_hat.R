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
omega_hat <- function(G, delta = NULL, eps = NULL, support = FALSE, sp = 0,
                      max_iter = 100, rel_tol = 1e-7, abs_tol = 1e-3, 
                      verbose = FALSE) {
  
  if (is.null(delta) && is.null(eps)) {
    # lambda <- .LambdaNR(t(G), 
    #                     max_iter = max_iter, 
    #                     rel_tol = rel_tol, 
    #                     support = support, 
    #                     verbose = verbose)
    # return(.OmegaHat(t(G), lambda = lambda, support = support))
    return(.OmegaHat(G = t(G), 
                     max_iter = max_iter,
                     rel_tol = rel_tol,
                     support = support,
                     verbose = verbose))
  }
  else {
    # Note: inital omegas obtained from non-censored optimization
    # lambda <- .LambdaNR(t(G), 
    #                     max_iter = max_iter, 
    #                     rel_tol = rel_tol, 
    #                     support = support, 
    #                     verbose = verbose)
    omegas_init <- .OmegaHat(G = t(G), 
                             max_iter = max_iter, 
                             rel_tol = rel_tol, 
                             support = support, 
                             verbose = verbose)
    
    if (any(is.nan(omegas_init))) {
      # stop("Initial omegas are nans.")
      return(rep(NaN, length(delta)))
    }
    else {
      if (sp == 0) {
        return(.OmegaHatEM(G = t(G), 
                           omegas_init = omegas_init, 
                           deltas = delta, 
                           epsilons = eps,
                           max_iter = max_iter, rel_tol = rel_tol, abs_tol = abs_tol,
                           support = support, 
                           verbose = verbose))
      }
      else {
        return(.OmegaHatEMSmooth(G = t(G), 
                                 omegas_init = omegas_init, 
                                 deltas = delta, 
                                 epsilons = eps,
                                 sp = sp,
                                 max_iter = max_iter, rel_tol = rel_tol, abs_tol = abs_tol,
                                 support = support, 
                                 erbose = verbose))
      }
    }
  }
}
