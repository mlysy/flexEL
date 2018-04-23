#' Evaluate the G matrix for mean regression (location-scale model)
#'
#' @param y Length-\code{nObs} vector of responses.
#' @param X \code{nObs x nBet} matrix of covariates.
#' @param Z \code{nObs x nGam} matrix of covariates.
#' @param beta Length-\code{nBet} vector of coefficients in location model.
#' @param gamma Length-\code{nGam} vector of coefficients in location model.
#' @param sig2 A scalar where sqrt(sig2) is the scale parameter for the error term.
#' @details Returns the G matrix for location-scale mean regression model. 
#' @return G matrix for location-scale mean regression model. 
#' @export mrls.evalG 
# mrls.evalG <- function(y, X, Z, beta, gamma) {
mrls.evalG <- function(y, X, Z, beta, gamma, sig2) {
  if (!is.vector(y)) stop("y should be a vector.") # TODO: allow y to be 1d matrix too
  if (nrow(X) != length(y)) stop("y and X have inconsistent dimensions.")
  if (nrow(Z) != length(y)) stop("y and Z have inconsistent dimensions.")
  # G <- .MeanRegLS_evalG(y,t(X),t(Z),beta,gamma)
  G <- .MeanRegLS_evalG(y,t(X),t(Z),beta,gamma,sig2)
  return(t(G))
}