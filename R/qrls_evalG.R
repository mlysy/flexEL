#' Evaluate the G matrix for a location-scale quantile regression model.
#'
#' @param y Response vector of length `n_obs`.
#' @param X Location covariate matrix of size `n_obs x n_bet`.
#' @param Z Scale covariate matrix of size `n_obs x n_gam`.
#' @param alphas Vector of `n_qts` quantile levels.
#' @param Beta Matrix of `n_bet x n_qts` location coefficients.
#' @param Gamma Matrix of `n_gam x n_qts` scale coefficients.
#' @param Nu Vector of `n_nu` initial values for the chain.
#' @example examples/qrls_evalG.R
#' @details ...
#' @return G matrix of size `for location-scale quantile regression model.
#' @export
qrls_evalG <- function(y, X, Z, alphas, Beta, Gamma, Sig2, Nu) {
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
  G <- .QuantRegLSEvalG(y, t(X), t(Z), alpha, Beta, Gamma, Sig2, Nu)
  return(t(G))
}
