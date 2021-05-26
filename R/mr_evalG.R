#' Evaluate the G matrix for a mean regression model.
#' 
#' @template args-y
#' @param X A numeric matrix of covariates of dimension \code{n_obs} x \code{n_bet} 
#'   where \code{n_obs} is the number of observations and \code{b_bet} is the number 
#'   of coefficients (length of \code{beta}).
#' @param beta A numeric vector of coefficients.
#' @details Assuming data were generated from 
#' ```
#' y_i = x_i'beta + eps_i, for i = 1, ..., n,
#' ```
#' where `eps_i`'s are ~iid `F(eps)`, with `E[eps] = 0` and `Var[eps] = 1`.
#' The estimating equations for `beta` can be derived based on the least squares method, 
#' that is,
#' ```
#' min_beta (y - X*beta)'(y - X*beta)
#' ```
#' The above minimization problem is equivalent to setting it's first derivative 
#' w.r.t `beta` to 0, i.e.,
#' ```
#' X'*(y - X*beta) = 0
#' ```
#' The left-hand-side of the above equation is the `G` matrix returned by this function.
#' @return A numeric matrix of dimension \code{n_obs} x \code{n_bet}.
#' @example examples/mr_evalG.R
#' @export mr_evalG
mr_evalG <- function(y, X, beta) {
  if (!is.vector(y)) stop("y should be a vector.") # TODO: allow y to be 1d matrix too
  if (nrow(X) != length(y)) stop("y and X have inconsistent dimensions.")
  G <- .MeanRegEvalG(y, t(X), beta)
  return(t(G))
}
