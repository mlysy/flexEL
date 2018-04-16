#' Evaluate the G matrix for quantile regression (location model)
#'
#' @param y Length-\code{nObs} vector of responses.
#' @param X \code{nObs x nEqs} matrix of covariates.
#' @param alphas a vector of quantile levels.
#' @param Beta \code{nEqs x nQts} matrix, each column is a vector of coefficients in location model.
#' @return G matrix for location quantile regression model. 
#' @details ...
#' @export qr.evalG
qr.evalG <- function(y, X, alphas, Beta) { 
  if (!is.vector(y)) stop("y should be a vector.") # TODO: allow y to be 1d matrix too
  if (nrow(X) != length(y)) stop("y and X have inconsistent dimensions.")
  # if input is for single quantile and Beta, Gamma are vectors, convert to matrix form
  if (is.vector(Beta)) Beta <- matrix(Beta,length(Beta),1)
  if (nrow(Beta) != ncol(X)) stop("X and Beta have inconsistent dimensions.")
  # add a first entry of alpha as the number of quantile levels
  alpha <- c(length(alphas), alphas) 
  G <- .QuantReg_evalG(y,t(X),alpha,Beta)
  return(t(G))
}