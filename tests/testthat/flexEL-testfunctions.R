#--- R version of C++ functions ------------------------------------------------

smooth_indicator <- function(eps1, eps2, s) {
  1/(1 + exp(s * (eps1 - eps2)))
}

expected_weights <- function(delta, epsilon, omega, smooth_s) {
  n_obs <- length(delta)
  n_obs2 <- length(omega)
  supp_adj <- n_obs2 == n_obs + 1
  weights <- rep(NA, n_obs2)
  if(supp_adj) {
    epsilon <- c(epsilon, -Inf)
    delta <- c(as.logical(delta), FALSE)
  }
  for(ii in 1:n_obs2) {
    omega_tilde <- smooth_indicator(epsilon[ii], epsilon, smooth_s)
    omega_tilde[ii] <- .5
    omega_tilde <- omega_tilde * omega
    omega_tilde <- omega_tilde / sum(omega_tilde)
    weights[ii] <- delta[ii] + sum((1-delta) * omega_tilde)
  }
  weights
}
