#' Evaluate the G matrix for the quantile regression model.
#'
#' @template args-y_X
#' @param alpha a vector of quantile levels.
#' @param Beta `n_eqs x n_qts` matrix, each column is a vector of coefficients in location model.
#' @example examples/qr_evalG.R
#' @return G matrix for location quantile regression model.
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
