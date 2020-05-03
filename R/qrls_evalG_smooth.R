#' Evaluate the G matrix for smoothed quantile regression location-scale model
#'
#' @param y Length-\code{n_obs} vector of responses.
#' @param X \code{n_obs x n_bet} matrix of covariates.
#' @param Z \code{n_obs x n_gam} matrix of covariates.
#' @param alphas a vector of quantile levels.
#' @param Beta \code{n_bet x n_qts} matrix, each column is a vector of coefficients in location function.
#' @param Gamma \code{n_gam x n_qts} matrix, each column is a vector of coefficients in scale function.
#' @param Nu A vector of quantile values of the same length as alpha.
#' @param s A positive scalar as smoothing parameter.
#' @example examples/qrls_evalG_smooth.R
#' @return G matrix for location-scale quantile regression model. 
#' @export qrls_evalG_smooth
qrls_evalG_smooth <- function(y, X, Z, alphas, Beta, Gamma, Sig2, Nu, sp = 10) { 
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
  G <- .QuantRegLSEvalGSmooth(y, t(X), t(Z), alpha, Beta, Gamma, Sig2, Nu, sp)
  return(t(G))
}