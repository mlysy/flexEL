#' Returns residuals of a location-scale model given the parameters.
#'
#' @param y Length-\code{nObs} vector of responses.
#' @param X Matrix of dimension \code{nObs x nBet} of location covariates.
#' @param Z Matrix of dimension \code{nObs x nGam} of scale covariates.
#' @param beta Length-\code{nObs} vector of responses.
#' @param rel_tol Relative tolerance of Newton-Raphson convergence.
#' @param verbose Display number of steps and tolerance criterion when algorithm terminates.
#' @return Length-\code{nObs} vector of residual epsilons.
#' @details ...
#' @export evalEpsilonsLS
evalEpsilonsLS <- function(y,X,Z,beta,gamma,sig2) {
  eps <- .evalEpsilonsLS(y,t(X),t(Z),beta,gamma,sig2)
  return(eps)
}