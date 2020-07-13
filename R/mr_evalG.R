#' Evaluate the G matrix for a mean regression model.
#'
#' @template args-y_X
#' @param beta Length-`n_eqs` vector of coefficients in location model.
#' @example examples/mr_evalG.R
#' @return G matrix for location mean regression model.
#' @export mr_evalG
mr_evalG <- function(y, X, beta) {
  if (!is.vector(y)) stop("y should be a vector.") # TODO: allow y to be 1d matrix too
  if (nrow(X) != length(y)) stop("y and X have inconsistent dimensions.")
  G <- .MeanRegEvalG(y,t(X),beta)
  return(t(G))
}
