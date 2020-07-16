#' Calculates log empirical likelihood.
#'
#' @template arg-G
#' @template arg-delta
#' @template arg-eps
#' @template arg-support
#' @example examples/logEL.R
#' @return A list containing the log EL value and the probability vector omega, or only the log EL value as a scalar.
#' @export logEL
logEL <- function(G, delta = NULL, eps = NULL, support = FALSE, 
                  max_iter = 100, rel_tol = 1e-7, abs_tol = 1e-3, 
                  return_omega = FALSE, verbose = FALSE) {
  
  if (is.null(delta) & is.null(eps)) {
    omega <- flexEL:::omega_hat(G = G, support = support, 
                                max_iter = max_iter, rel_tol = rel_tol, abs_tol = abs_tol, verbose = verbose)
  }
  else if (!is.null(delta) & !is.null(eps)) {
    omega <- flexEL:::omega_hat(G = G, delta = delta, eps = eps, support = support,
                                max_iter = max_iter, rel_tol = rel_tol, abs_tol = abs_tol, verbose = verbose)
  }
  else if (is.null(delta) + is.null(eps) == 1) {
    stop("`delta` and `eps` must be given together if one is given.")
  }
  
  if (return_omega) {
    return(list(log_el = .LogEL(omegas = omega, support = support),
                omega = omega))
  }
  else{
    return(.LogEL(omegas = omega, support = support))
  }
}
