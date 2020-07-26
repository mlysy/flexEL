#' Evaluate the G matrix for a quantile regression model.
#' 
#' @template args-y_X
#' @param alpha a vector of quantile levels.
#' @param Beta `n_eqs x n_qts` matrix, each column is a vector of coefficients in location model.
#' @details 
#' Assuming data were generated from 
#' ```
#' y_i = x_i'beta + eps_i, for i = 1, ..., n,
#' ```
#' where `eps_i`'s are ~iid `eps`, with `E[eps] = 0` and `Var[eps] = 1`. 
#' Quantile regression estimates the alpha-level quantile of the response variable, i.e., 
#' ```
#' Q_alpha(y | x_i) = x_i'beta, for i = 1, ..., n.
#' ```
#' where the 1st element of `x_i` is 1 and such that the 1st element of `beta` corresponds 
#' to the `alpha`-level quantile of `eps`.
#' @return G matrix for location quantile regression model.
#' @example examples/qr_evalG.R
#' @export qr_evalG
qr_evalG <- function(y, X, alpha, Beta) {
  
  if (!is.vector(y)) stop("y should be a vector.") # TODO: allow y to be 1d matrix too
  if (nrow(X) != length(y)) stop("y and X have inconsistent dimensions.")
  # if input is for single quantile and Beta, Gamma are vectors, convert to matrix form
  if (length(alpha) == 1 && is.vector(Beta)) Beta <- matrix(Beta, length(Beta), 1)
  if (length(alpha) > 1 && is.vector(Beta)) {
    stop("Parameters must be in matrix form when alpha has more than one entry.")
  }
  if (nrow(Beta) != ncol(X)) stop("X and Beta have inconsistent dimensions.")
  
  # the first entry of alpha passed to the C++ function is the number of quantile levels
  alpha <- c(length(alpha), alpha)
  G <- .QuantRegEvalG(y, t(X), alpha, Beta)
  return(t(G))
}
