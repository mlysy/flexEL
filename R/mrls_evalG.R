#' Evaluate the G matrix for a location-scale mean regression model.
#'
#' @template args-y_X
#' @template arg-Z
#' @param beta Length-`n_bet` vector of coefficients in location model.
#' @param gamma Length-`n_gam` vector of coefficients in location model.
#' @param sig2 A positive scalar whose square root is the scale parameter for the error term.
#' @example examples/mrls_evalG.R
#' @return G matrix for location-scale model.
#' @export mrls_evalG
mrls_evalG <- function(y, X, Z, beta, gamma, sig2) {
  if (!is.vector(y)) stop("y should be a vector.") # TODO: allow y to be 1d matrix too
  if (nrow(X) != length(y)) stop("y and X have inconsistent dimensions.")
  if (nrow(Z) != length(y)) stop("y and Z have inconsistent dimensions.")
  # G <- .MeanRegLS_evalG(y,t(X),t(Z),beta,gamma)
  G <- .MeanRegLSEvalG(y,t(X),t(Z),beta,gamma,sig2)
  return(t(G))
}
