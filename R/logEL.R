#' Calculates log empirical likelihood.
#'
#' @template arg-G
#' @template arg-delta
#' @template arg-eps
#' @template arg-support
#' @template arg-sp
#' @example examples/logEL.R
#' @return A list containing the log EL value and the probability vector omega, or only the log EL value as a scalar.
#' @export logEL
logEL <- function(G, delta = NULL, eps = NULL, support = FALSE, sp = 0, 
                  max_iter = 100, rel_tol = 1e-7, abs_tol = 1e-3, 
                  return_omega = FALSE, verbose = FALSE) {
  
  if (is.null(delta) & is.null(eps)) {
    omega <- flexEL:::omega_hat(G = G, support = support, 
                                max_iter = max_iter, rel_tol = rel_tol, abs_tol = abs_tol, verbose = verbose)
    if (return_omega) {
      return(list(log_el = .LogEL(omegas = omega, support = support),
                  omega = omega))
    }
    else{
      return(.LogEL(omegas = omega, support = support))
    }
  }
  else if (!is.null(delta) & !is.null(eps)) {
    omega <- flexEL:::omega_hat(G = G, delta = delta, eps = eps, support = support,
                                max_iter = max_iter, rel_tol = rel_tol, abs_tol = abs_tol, verbose = verbose)
    if (sp == 0) {
      if (return_omega) {
        return(list(log_el = .LogELCens(omegas = omega, epsilons = eps, deltas = delta, support = support),
                    omega = omega))
      }
      else{
        return(.LogELCens(omegas = omega, epsilons = eps, deltas = delta, support = support))
      }
    }
    else if (sp > 0) {
      if (return_omega) {
        return(list(log_el = .LogELSmooth(omegas = omega, epsilons = eps, deltas = delta, 
                                          sp = sp, support = support),
                    omega = omega))
      }
      else{
        return(.LogELSmooth(omegas = omega, epsilons = eps, deltas = delta, 
                            sp = sp, support = support))
      }
    }
    else {
      stop("`sp` must be a positive number.")
    }
  }
  else if (is.null(delta) + is.null(eps) == 1) {
    stop("`delta` and `eps` must be given together if one is given.")
  }
}
