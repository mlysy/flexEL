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
  
  omega <- flexEL:::omega_hat(G, delta, eps, support,
                              max_iter, rel_tol, abs_tol, verbose)
  if (return_omega) {
    return(list(log_el = .LogEL(omega, support),
                omega = omega))
  }
  else{
    return(.LogEL(omega, support))
  }
}
