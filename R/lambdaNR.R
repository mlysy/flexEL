#' Solve the inner loop optimization of an EL function.
#'
#' @param G \code{nObs x nEqns} matrix of constraints.
#' @param lambda0 Length-\code{nEq} vector of starting values.
#' @param max_iter Maximum number of Newton-Raphson steps.
#' @param eps Relative tolerance of Newton-Raphson convergence.
#' @param verbose Display number of steps and tolerance criterion when algorithm terminates.
#' @return Length-\code{nEq} vector corresponding to the solution of the optimization problem.
#' @details The inner-loop optimization of EL is ...
#' @export
lambdaNR <- function(G, lambda0, nObs, nEqs,
                     max_iter = 100, eps = 1e-7, verbose = FALSE) {
  # input checks
  nEqs <- ncol(G)
  nObs <- nrow(G)
  if(length(lambda0) != nEqs) {
    stop("G and lambda0 have incompatible dimensions.")
  }
  ans <- .lambdaNR(G = t(G),
                   maxIter = max_iter, eps = eps, verbose = verbose)
  if(ans$convergence) {
    ans <- ans$lambda
  } else {
    ans <- rep(NA, nEqs)
  }
  ans
}
