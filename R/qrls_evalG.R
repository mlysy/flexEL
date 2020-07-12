#' Evaluate the G matrix for a location-scale quantile regression model.
#'
#' @template args-y
#' @template args-X
#' @template args-Z
#' @param alphas A vector of quantile levels.
#' @param Beta An `n_bet x n_qts` matrix, each column is a vector of coefficients in location function.
#' @param Gamma An `n_gam x n_qts` matrix, each column is a vector of coefficients in scale function.
#' @param Nu A vector of quantile values of the same length as alpha.
#' @param sp A positive scalar as smoothing parameter.
#' @example examples/qrls_evalG.R
#' @return G matrix of size `for location-scale quantile regression model.
#' @export qrls_evalG
qrls_evalG <- function(y, X, Z, alphas, Beta, Gamma, Sig2, Nu, sp = 0) {
  
  # checks
  if (sp <= 0) stop("s must be positive.")
  if (!is.vector(y)) stop("y should be a vector.") # TODO: allow y to be 1d matrix too
  if (nrow(X) != length(y)) stop("y and X have inconsistent dimensions.")
  if (nrow(Z) != length(y)) stop("y and Z have inconsistent dimensions.")
  
  # if input is for single quantile and Beta, Gamma are vectors, convert to matrix form
  if (length(alphas) == 1 && is.vector(Beta)) Beta <- matrix(Beta,length(Beta),1)
  if (length(alphas) == 1 && is.vector(Gamma)) Gamma <- matrix(Gamma,length(Gamma),1)
  if (length(alphas) > 1 && (is.vector(Beta) || is.vector(Gamma))) {
    stop("Parameters must be in matrix form when alphas has more than one entry.")
  }
  
  # the first entry of alpha passed to the C++ function is the number of quantile levels
  alpha <- c(length(alphas), alphas)
  if (sp == 0) {
    return(.QuantRegLSEvalG(y, t(X), t(Z), alpha, Beta, Gamma, Sig2, Nu))
  }
  else {
    return(QuantRegLSEvalGSmooth(y, t(X), t(Z), alpha, Beta, Gamma, Sig2, Nu, sp))
  }
}
