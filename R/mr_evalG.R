#' Evaluate the G matrix for mean regression (location model)
#'
#' @param y Length-\code{n_obs} vector of responses.
#' @param X \code{n_obs x n_eqs} matrix of covariates.
#' @param beta Length-\code{n_eqs} vector of coefficients in location model.
#' @return G matrix for location mean regression model. 
#' @details ...
#' @export mr_evalG
mr_evalG <- function(y, X, beta) { 
  if (!is.vector(y)) stop("y should be a vector.") # TODO: allow y to be 1d matrix too
  if (nrow(X) != length(y)) stop("y and X have inconsistent dimensions.")
  G <- .MeanRegEvalG(y,t(X),beta)
  return(t(G))
}