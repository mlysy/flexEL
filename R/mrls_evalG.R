#' Evaluate the G matrix for a location-scale mean regression model.
#'
#' @template args-y_X_Z
#' @param beta A numeric vector of coefficients in the location function.
#' @param gamma A numeric vector of coefficients in the scale function.
#' @param sig2 A positive scalar whose square root is the scale parameter in the scale function.
#' @details Assuming data were generated from 
#' ```
#' y_i = x_i'beta + sqrt(sig2) * exp(z_i'gamma) * eps_i, for i = 1, ..., n,
#' ```
#' where `eps_i`'s are ~iid `F(eps)`, with `E[eps] = 0` and `Var[eps] = 1`.
#' In this location-scale model, `x_i'beta` is the location function and `sqrt(sig2) * exp(z_i'gamma)` 
#' is the scale function. The `G` matrix is calculated by stacking together the first derivative of the 
#' quasi-likelihood function w.r.t `beta`, `gamma`, and `sig2`. 
#' More details of this setup can be found in the package vignette: 
#' \code{vignette("help", package = "flexEL")}.
#' @return A numeric matrix of dimension \code{n_obs} x \code{n_bet + n_gam + 1}.
#' @example examples/mrls_evalG.R
#' @export mrls_evalG
mrls_evalG <- function(y, X, Z, beta, gamma, sig2) {
  if (!is.vector(y)) stop("y should be a vector.") # TODO: allow y to be 1d matrix too
  if (nrow(X) != length(y)) stop("y and X have inconsistent dimensions.")
  if (nrow(Z) != length(y)) stop("y and Z have inconsistent dimensions.")
  G <- MeanRegLS_evalG(y, t(X), t(Z), beta, gamma, sig2)
  return(t(G))
}
