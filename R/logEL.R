#' Calculates log empirical likelihood.
#'
#' @template arg-G
#' @template arg-delta
#' @template arg-eps
#' @template arg-support
#' @template arg-sp
#' @template args-lambda_precision
#' @template arg-abs_tol
#' @param return_omega A boolean indicating whether to return the probability vector omega used to evaluate log EL.
#' @param return_dldG A boolean indicating whether to return the gradient matrix dldG of log EL w.r.t. G.
#' @template arg-verbose
#' @example examples/logEL.R
#' @return A list containing the log EL value, the probability vector omega, and a matrix 
#'   containing the gradient of log EL function with respect to G matrix evaluated at G, 
#'   or only the log EL value as a scalar.
#' @export logEL
logEL <- function(G, lambda0 = rep(0, ncol(G)),
                  delta = NULL, eps = NULL, support = FALSE, sp = 0, 
                  max_iter = 100, rel_tol = 1e-7, abs_tol = 1e-3, 
                  return_omega = FALSE, return_dldG = FALSE, verbose = FALSE) {
  
  if (is.null(delta) & is.null(eps)) { # no censor
    
    if (return_dldG) {
      return(.LogELGrad(G = t(G), 
                        lambda0 = lambda0,
                        max_iter = max_iter, 
                        rel_tol = rel_tol, 
                        support = support, 
                        verbose = verbose))
    }
    else {
      omega <- omega_hat(G = G, 
                         lambda0 = lambda0,
                         support = support, 
                         max_iter = max_iter, 
                         rel_tol = rel_tol, 
                         abs_tol = abs_tol, 
                         verbose = verbose)
      if (return_omega) {
        return(list(log_el = .LogEL(omegas = omega, support = support),
                    omega = omega))
      }
      else {
        return(.LogEL(omegas = omega, support = support))
      }
    }
  }
  else if (!is.null(delta) & !is.null(eps)) { # right-censor
    
    if (return_dldG) {
      stop("Currently does not support gradient calculation for right-censored data.")
    }
    
    omega <- omega_hat(G = G, lambda0 = lambda0,
                       delta = delta, eps = eps, support = support,
                       max_iter = max_iter, rel_tol = rel_tol, abs_tol = abs_tol, 
                       verbose = verbose)
    if (sp == 0) {
      if (return_omega) {
        return(list(log_el = .LogELCens(omegas = omega, epsilons = eps, deltas = delta, 
                                        support = support),
                    omega = omega))
      }
      else{
        return(.LogELCens(omegas = omega, deltas = delta, epsilons = eps, support = support))
      }
    }
    else if (sp > 0) {
      if (return_omega) {
        return(list(log_el = .LogELSmooth(omegas = omega, deltas = delta, epsilons = eps, 
                                          sp = sp, support = support),
                    omega = omega))
      }
      else{
        return(.LogELSmooth(omegas = omega, deltas = delta, epsilons = eps, 
                            sp = sp, support = support))
      }
    }
    else {
      stop("`sp` must be a positive number.")
    }
  }
  else if (is.null(delta) + is.null(eps) == 1) {
    stop("`delta` and `eps` must both provided if one is provided")
  }
}
