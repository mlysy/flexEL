#' Evaluate the G matrix for a location-scale quantile regression model.
#'
#' @template args-y_X_Z
#' @param alpha A numeric vector of quantile levels of length \code{n_qts}.
#' @param Beta A numeric matrix of dimension \code{n_bet} x \code{n_qts}. Each column of \code{Beta} 
#'   is a vector of coefficients in the location function.
#' @param Gamma A An numeric matrix of dimension \code{n_gam} x \code{n_qts}. Each column of \code{Beta} 
#'   is a vector of coefficients in the scale function.
#' @param Sig2 A positive numeric vector of length \code{n_qts} where each element' square root is 
#'   the scale parameter in the scale function corresponding to the quantile level.
#' @param Nu A numeric vector of quantile values of length \code{n_qts} corresponding to each alpha.
#' @param sp A positive scalar as the smoothing parameter. If `sp = 0`, then no smoothing is performed.
#' @details Assuming data were generated from 
#' ```
#' y_i = x_i'beta + sigma * exp(z_i'gamma) * eps_i, for i = 1, ..., n,
#' ```
#' where `eps_i`'s are ~iid `F(eps)`, with `E[eps] = 0` and `Var[eps] = 1`. 
#' Quantile regression estimates the alpha-level quantile of the response variable, i.e., 
#' ```
#' Q_alpha(y | x_i, z_i) = x_i'beta + sigma * exp(z_i'gamma) * nu_alpha,
#' for i = 1, ..., n
#' ```
#' where `nu_alpha` is the alpha-level quantile value of `eps`. 
#' In this location-scale model, `x_i'beta` is the location function and `sqrt(sig2) * exp(z_i'gamma)` 
#' is the scale function. The `G` matrix is calculated by stacking together the first derivative of the 
#' quasi-likelihood function w.r.t `beta`, `gamma`, `sig2`, and the first derivative w.r.t `nu_alpha` 
#' of the check function introduced by Koenker and Bassett (1978),
#' ```
#' rho_alpha(u) = u * (alpha - 1{u <= 0})
#' ```
#' where `alpha` is the quantile level and `1{}` is the indicator function which returns 1 if the
#' condition is true and 0 otherwise.
#' More details of this setup can be found in the package vignette: 
#' \code{vignette("help", package = "flexEL")}.
#' @references G. Basset and R. Koenker. Regression quantiles. Econometrica, 46(1):33–50, 1978.
#' @return A numeric matrix of dimension \code{n_obs} x \code{n_bet + n_gam + 2}.
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
