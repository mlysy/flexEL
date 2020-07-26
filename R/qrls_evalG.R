#' Evaluate the G matrix for a location-scale quantile regression model.
#'
#' @template args-y_X
#' @template arg-Z
#' @param alpha A length-`n_qts` numeric vector of quantile levels.
#' @param Beta An `n_bet x n_qts` matrix, each column is a vector of coefficients in location function.
#' @param Gamma An `n_gam x n_qts` matrix, each column is a vector of coefficients in scale function.
#' @param Sig2 A positive scalar whose square root is the scale parameter for the error term.
#' @param Nu A length-`n_qts` numeric vector of quantile values corresponding to each alpha.
#' @param sp A positive scalar as smoothing parameter.
#' @details Assuming data were generated from 
#' ```
#' y_i = x_i'beta + sigma * exp(z_i'gamma) * eps_i, for i = 1, ..., n,
#' ```
#' where `eps_i`'s are ~iid `eps`, with `E[eps] = 0` and `Var[eps] = 1`. 
#' Quantile regression estimates the alpha-level quantile of the response variable, i.e., 
#' ```
#' Q_alpha(y | x_i, z_i) = x_i'beta + sigma * exp(z_i'gamma) * nu_alpha,
#' for i = 1, ..., n.
#' ```
#' where `nu_alpha` is the alpha-level quantile value of `eps`. 
#' Neither `x_i` nor `z_i` should have a constant term.
#' @return G matrix of size `for location-scale quantile regression model.
#' @example examples/qrls_evalG.R
#' @export qrls_evalG
qrls_evalG <- function(y, X, Z, alpha, Beta, Gamma, Sig2, Nu, sp = 0) {
  
  # checks
  if (sp < 0) stop("s must be positive.")
  if (!is.vector(y)) stop("y should be a vector.") # TODO: allow y to be 1d matrix too
  if (nrow(X) != length(y)) stop("y and X have inconsistent dimensions.")
  if (nrow(Z) != length(y)) stop("y and Z have inconsistent dimensions.")
  
  # if input is for single quantile and Beta, Gamma are vectors, convert to matrix form
  if (length(alpha) == 1 && is.vector(Beta)) Beta <- matrix(Beta,length(Beta),1)
  if (length(alpha) == 1 && is.vector(Gamma)) Gamma <- matrix(Gamma,length(Gamma),1)
  if (length(alpha) > 1 && (is.vector(Beta) || is.vector(Gamma))) {
    stop("Parameters must be in matrix form when alpha has more than one entry.")
  }
  
  # the first entry of alpha passed to the C++ function is the number of quantile levels
  alpha <- c(length(alpha), alpha)
  if (sp == 0) {
    return(t(.QuantRegLSEvalG(y, t(X), t(Z), alpha, Beta, Gamma, Sig2, Nu)))
  }
  else {
    return(t(.QuantRegLSEvalGSmooth(y, t(X), t(Z), alpha, Beta, Gamma, Sig2, Nu, sp)))
  }
}
