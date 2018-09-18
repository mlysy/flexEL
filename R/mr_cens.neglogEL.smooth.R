#' Calculate the negative log EL with smoothing for mean regression under censoring
#'
#' @template args-y
#' @template args-X
#' @template args-delta
#' @param beta Length-\code{nBet} vector of coefficients in mean regression model.
#' @template args-sp
#' @return Scalar of negative smoothed log EL.
#' @details ...
#' @export mr_cens.neglogEL.smooth
mr_cens.neglogEL.smooth <- function(y, X, delta, beta, sp=10) {
  G <- mr.evalG(y, X, beta)
  eps <- c(y-X%*%beta)
  omegas <- omega.hat.EM.smooth(G,delta,eps,sp)
  if (!anyNA(omegas)) {
    res <- -logEL.smooth(omegas,eps,delta)
    # gradlist <- mr.deltaG_R(y, X, beta)
    # grad <- -logELCensgrad_R(omegas, deltas, epsilons, lambda, gradlist, weights) # negative gradient
    # attr(res, "gradient") <- grad
    # attr(res, "gradient") <- grad(-logEL.smooth_R,)
  }
  else {
    # TODO: if not converged, what should be the gradient...??
    res <- Inf
    attr(res, "gradient") <- Inf
  }
  return(res)
}