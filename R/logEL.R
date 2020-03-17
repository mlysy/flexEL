#' log empirical likelihood
#'
#' @template args-omega
#' @template args-eps
#' @template args-delta
#' @return Log empirical likelihood of the input beta
#' @details ...
#' @export logEL
logEL <- function(omega, eps, delta, support = FALSE) {
  if (missing(eps) && missing(delta)) {
    logel <- .LogEL(omega, support)
  }
  else {
    logel <- .LogELC(omega,eps,delta,support)
  }
  return(logel)
}

# logEL <- function(G, deltas, epsilons, 
#                   max_iter = 100, rel_tol = 1e-7, verbose = FALSE) {
#   # non-censoing case
#   if (missing(deltas) && missing(epsilons)) {
#     lambda <- .lambdaNR(t(G), maxIter = max_iter, relTol = rel_tol, verbose = verbose)
#     omegas <- .omega.hat(t(G), lambda = lambda)
#     logel <- .logEL(omegas)
#     # .logEL(G = t(G), maxIter = max_iter, relTol = rel_tol, verbose = verbose)
#   }
#   # censoring case
#   else {
#     # input checks
#     if(nrow(G) != length(deltas)) {
#       stop("G and deltas have inconsistent dimensions.")
#     }
#     if(nrow(G) != length(epsilons)) {
#       stop("G and epsilons have inconsistent dimensions.")
#     }
#     if(length(deltas) != length(epsilons)) {
#       stop("deltas and epsilons have inconsistent lengths.")
#     }
#     # Note: inital omegas obtained from non-censored optimization
#     lambda <- .lambdaNR(t(G), maxIter = max_iter, relTol = rel_tol, verbose = verbose)
#     omegasInit <- .omega.hat(t(G), lambda)
#     omegas <- .omega.hat.EM(omegasInit, t(G), deltas, epsilons, max_iter, rel_tol, verbose)
#     # if (any(is.nan(omegas))) return(-Inf)
#     logel <-.logELC(omegas = omegas, G = t(G), deltas = deltas, epsilons = epsilons)
#   }
#   return(logel)
# }
