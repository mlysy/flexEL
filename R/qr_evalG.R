#' Evaluate the G matrix for a quantile regression model.
#' 
#' @template args-y
#' @param X A numeric matrix of covariates of dimension \code{n_obs} x \code{n_bet} 
#'   where \code{n_obs} is the number of observations and \code{n_bet} is the number 
#'   of coefficients (length of \code{beta}). `X` should contain a column of `1`s and 
#'   the corresponding element of `beta` (a column of `Beta`) corresponds to the 
#'   quantile value of `eps` at that quantile level.
#' @param alpha a vector of quantile levels.
#' @param Beta A numeric matrix of dimension \code{n_bet} x \code{n_qts}, each column is a 
#'   vector of coefficients in location model corresponding to the quantile level. 
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
#' where the 1st element of `x_i` should be 1 and  the 1st element of `beta` corresponds 
#' to the `alpha`-level quantile of `eps`.
#' The estimating equation can be derived based on the minimizing "check function" introduced by 
#' Koenker and Bassett (1978),
#' ```
#' rho_alpha(u) = u * (alpha - 1{u <= 0})
#' ```
#' where `alpha` is the quantile level, `1{}` is the indicator function which returns 1 if the
#' condition is true and 0 otherwise, and `u_i = y_i - x_i'beta`. That is,
#' ```
#' min_beta rho(y - X*beta)
#' ```
#' The above minimization problem is equavalent to setting it's first derivative 
#' w.r.t `beta` to 0, i.e.,
#' ```
#' phi(y - X*beta) = 0
#' ```
#' where `phi` is the first derivative of `rho`.
#' The left-hand-side of the above equation is the `G` matrix returned by this function.
#' @references G. Basset and R. Koenker. Regression quantiles. Econometrica, 46(1):33–50, 1978.
#' @return A numeric matrix of dimension \code{n_obs} x \code{n_bet}.
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
