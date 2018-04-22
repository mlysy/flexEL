#' Evaluate the G matrix for mean regression (location-scale model)
#'
#' @param y Length-\code{nObs} vector of responses.
#' @param X \code{nObs x nBet} matrix of covariates.
#' @param Z \code{nObs x nGam} matrix of covariates.
#' @param alphas a vector of quantile levels.
#' @param Beta \code{nBet x nQts} matrix, each column is a vector of coefficients in location function.
#' @param Gamma \code{nGam x nQts} matrix, each column is a vector of coefficients in scale function.
#' @param Nu Length-\code{numNu} vector of initial value for the chain.
#' @return G matrix for location-scale quantile regression model. 
#' @export qrls.evalG
qrls.evalG <- function(y, X, Z, alphas, Beta, Gamma, Nu) { 
  if (!is.vector(y)) stop("y should be a vector.") # TODO: allow y to be 1d matrix too
  if (nrow(X) != length(y)) stop("y and X have inconsistent dimensions.")
  if (nrow(Z) != length(y)) stop("y and Z have inconsistent dimensions.")
  # if input is for single quantile and Beta, Gamma are vectors, convert to matrix form
  if (length(alphas) == 1 && is.vector(Beta)) Beta <- matrix(Beta,length(Beta),1)
  if (length(alphas) == 1 && is.vector(Gamma)) Gamma <- matrix(Gamma,length(Gamma),1)
  if (length(alphas) > 1 && (is.vector(Beta) || is.vector(Gamma))) {
    stop("Parameters must be in matrix form when alphas has more than one entry.")
  }
  alpha <- c(length(alphas), alphas) 
  G <- .QuantRegLS_evalG(y, t(X), t(Z), alpha, Beta, Gamma, Nu)
  return(t(G))
}