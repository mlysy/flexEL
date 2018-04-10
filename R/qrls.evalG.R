#' Evaluate the G matrix for mean regression (location-scale model)
#'
#' @param y Length-\code{nObs} vector of responses.
#' @param X \code{nObs x nBet} matrix of covariates.
#' @param Z \code{nObs x nGam} matrix of covariates.
#' @param alpha a scalar of quantile leve.
#' @param beta Length-\code{nBet} vector of coefficients in location model.
#' @param gamma Length-\code{nGam} vector of coefficients in location model.
#' @details Returns the G matrix for location-scale quantile regression model. 
#' @export
qrls.evalG <- function(y, X, Z, alpha, beta, gamma) { 
  if (!is.vector(y)) stop("y should be a vector.") # TODO: allow y to be 1d matrix too
  if (nrow(X) != length(y)) stop("y and X have inconsistent dimensions.")
  if (nrow(Z) != length(y)) stop("y and Z have inconsistent dimensions.")
  G <- .QuantRegLS_evalG(y, t(X), t(Z), alpha, beta, gamma)
  return(t(G))
}