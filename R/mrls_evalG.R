#' Evaluate the G matrix for a location-scale mean regression model.
#'
#' @template args-y_X
#' @template arg-Z
#' @param beta Length-`n_bet` vector of coefficients in location model.
#' @param gamma Length-`n_gam` vector of coefficients in location model.
#' @param sig2 A positive scalar whose square root is the scale parameter for the error term.
#' @details Assuming data were generated from 
#' ```
#' y_i = x_i'beta + sigma * exp(z_i'gamma) * eps_i, for i = 1, ..., n,
#' ```
#' where `eps_i`'s are ~iid `eps`, with `E[eps] = 0` and `Var[eps] = 1`.
#' @return G matrix for location-scale model.
#' @example examples/mrls_evalG.R
#' @export mrls_evalG
mrls_evalG <- function(y, X, Z, beta, gamma, sig2) {
  if (!is.vector(y)) stop("y should be a vector.") # TODO: allow y to be 1d matrix too
  if (nrow(X) != length(y)) stop("y and X have inconsistent dimensions.")
  if (nrow(Z) != length(y)) stop("y and Z have inconsistent dimensions.")
  # G <- .MeanRegLS_evalG(y,t(X),t(Z),beta,gamma)
  G <- .MeanRegLSEvalG(y,t(X),t(Z),beta,gamma,sig2)
  return(t(G))
}
