#' Evaluate the G matrix for a mean regression model.
#' 
#' @template args-y_X
#' @param beta Length-`n_eqs` vector of coefficients in location model.
#' @details Assuming data were generated from 
#' ```
#' y_i = x_i'beta + eps_i, for i = 1, ..., n,
#' ```
#' where `eps_i`'s are ~iid `eps`, with `E[eps] = 0` and `Var[eps] = 1`.
#' @return G matrix for location mean regression model.
#' @example examples/mr_evalG.R
#' @export mr_evalG
mr_evalG <- function(y, X, beta) {
  if (!is.vector(y)) stop("y should be a vector.") # TODO: allow y to be 1d matrix too
  if (nrow(X) != length(y)) stop("y and X have inconsistent dimensions.")
  G <- .MeanRegEvalG(y,t(X),beta)
  return(t(G))
}
