/// @file logistic_el.stan
///
/// @brief Logistic regression by empirical likelihood 
///       using Bartlett moment conditions.

functions {
#include "include/stanEL.stan"
}

data {
  int<lower=1> n_obs; // Number of observations.
  int<lower=1> n_cov; // Number of covariates.
  vector[n_obs] y; // Number of successes per area/observation.
  vector[n_obs] offset; // Offset term.
  matrix[n_obs, n_cov] X; // Common covariate matrix.
  real<lower=0> beta_sd; // Common standard deviation for `beta`.
  // EL options:
  int<lower=1> max_iter;
  real<lower=0> rel_tol;
  int<lower=0,upper=1> supp_adj;
}

parameters {
  vector[n_cov] beta_raw;
}
transformed parameters {
  vector[n_cov] beta = beta_sd * beta_raw;
}
model {
  int dnc_inf = 1;
  // G matrix:
  vector[n_obs] Xb = X * beta + offset;
  vector[n_obs] rho = inv_logit(Xb);
  matrix[2, n_obs] G = rep_matrix(0, 2, n_obs);
  G[1] = (y - rho)';
  G[2] = ((y - rho) .* (y - rho) ./ (rho .* (1 - rho)) - 1)';
  target += logEL_opt(G, max_iter, rel_tol, supp_adj, dnc_inf);
  // prior:
  beta_raw ~ std_normal();
}



