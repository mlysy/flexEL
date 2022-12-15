
#' Calculate derivative of G w.r.t theta for location-scale mean regression model
#' 
#' @return A list of length `n_obs`.
#' 
#' @export
mrls_dGdt <- function(y, X, Z, beta, gamma, sig2) {
  n_obs <- length(y)
  n_eqs <- length(beta) + length(gamma) + 1
  dGdt <- MeanRegLS_dGdt(y, t(X), t(Z), beta, gamma, sig2)
  dGdt_lst <- split(dGdt, rep(seq(from = 1, to = n_obs), each = n_eqs))
  lapply(dGdt_lst, function(ll) {
    matrix(ll, nrow = n_eqs)
  })
}