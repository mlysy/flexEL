#' Calculate derivative of G w.r.t theta for mean regression model
#' 
#' @template args-y
#' @param X A numeric matrix of covariates of dimension \code{n_obs} x \code{n_bet} 
#'   where \code{n_obs} is the number of observations and \code{b_bet} is the number 
#'   of coefficients (length of \code{beta}).
#' @param beta A numeric vector of coefficients.
#' 
#' @return A list of length `n_obs`.
#' 
#' @export
mr_dGdt <- function(y, X, beta) {
  n_obs <- length(y)
  n_eqs <- length(beta)
  dGdt <- MeanReg_dGdt(y, t(X), beta)
  dGdt_lst <- split(dGdt, rep(seq(from = 1, to = n_obs), each = n_eqs))
  lapply(dGdt_lst, function(ll) {
    matrix(ll, nrow = n_eqs)
  })
}