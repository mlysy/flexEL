
source("./tests/testthat/el_rfuns.R")

# n <- 500 # number of observations
# b <- c(0.5,1) # beta_0 = 0.5, beta_1 = 1
# # eps <- gen_nct_eps(n, df = 20, ncp = 1) # a re-centered right-skewed non-central t distribution
# eps <- rnorm(n)
# X <- cbind(1, rnorm(n)) # n x 2 covariate matrix (intercept included)
# y <- X %*% b + eps

# # trying out for support adjust
# mr_neglogEL_adj_R <- function(y, X, beta) {
#   n <- length(y)
#   p <- length(beta)
#   G <- mr_evalG_R(y, X, beta)
#   neglogel <- - flexEL::logEL(G, support = TRUE)
#   # G <- adjG_R(G)
#   lambda <- flexEL::lambdaNR(G, support = TRUE)
#   omega <- flexEL:::omega_hat(G, support = TRUE)
#   # message(omega[n+1])
#   dldGadj <- logEL_dldG_R(lambda, omega[1:n]) + omega[n+1]*0.5*log(n)/n * rep(1,ncol = n)
#   dGdb <- mr_dGdt_R(y, X, beta)
#   grad_mat <- matrix(NA, nrow = nrow(dldGadj), ncol = ncol(dldGadj))
#   for (ii in 1:ncol(dldGadj)) {
#     grad_mat[,ii] <- dGdb[[ii]] %*% dldGadj[,ii]
#   }
#   grad <- rowSums(grad_mat)
#   attr(neglogel, "gradient") <- grad
#   return(neglogel)
# }

mr_neglogEL_R <- function(G) {
  lambda <- flexEL::lambdaNR(G = G, rel_tol = 1e-4, support = FALSE)
  # omega <- c(1/(1-t(lambda) %*% t(G)) / sum(1/(1-t(lambda) %*% t(G))))
  omega <- 1/nrow(G) * 1/(1 - c(G %*% lambda))
  neglogel <- -sum(log(omega))
  dldG <- logEL_dldG_R(lambda, omega)
  attr(neglogel, "gradient") <- dldG
  return(neglogel)
}

# coef(lm(y ~ X-1))
# bb <- c(0.5, 1.2)
# G <- mr_evalG_R(y, X, bb)

n_obs <- sample(10:30, 1)
n_eq <- sample(2:4, 1)
G <- matrix(rnorm(n_obs*n_eq), n_obs, n_eq)

res <- mr_neglogEL_R(G)
dldG_re <- attr(res, "gradient")
dldG_nd <- matrix(numDeriv::grad(mr_neglogEL_R, G), nrow = nrow(G), ncol = ncol(G))
head(dldG_re)
head(dldG_nd)

mr_neglogEL_adj_R <- function(G) {
  n <- nrow(G)
  lambda <- flexEL::lambdaNR(G = G, rel_tol = 1e-4, support = TRUE)
  omega <- flexEL:::omega_hat(G = G, support = TRUE)
  neglogel <- -sum(log(omega))
  an <- max(1, 0.5*log(n))
  dldGadj <- logEL_dldG_R(lambda, omega[1:n])/n*(n+1) + (n+1) * omega[n+1]*an/n * rep(1, n) %*% t(lambda)
  attr(neglogel, "gradient") <- dldGadj
  return(neglogel)
}

# n_obs <- 10
# n_eq <- 2
G <- matrix(rnorm(n_obs*n_eq), n_obs, n_eq)

# G <- mr_evalG_R(y, X, b)
# res <- mr_neglogEL_adj_R(G)

res <- mr_neglogEL_adj_R(G)
dldG_re <- attr(res, "gradient")
dldG_nd <- matrix(numDeriv::grad(mr_neglogEL_adj_R, G), nrow = nrow(G), ncol = ncol(G))
head(dldG_re)
head(dldG_nd)

#--- simple test following qin-lawless 1994 ------------------------------------

require(flexEL)
require(numDeriv)

n_obs <- 10
n_eq <- 3
G <- matrix(rnorm(n_obs*n_eq), n_obs, n_eq)

lambda <- lambdaNR(G, verbose = TRUE)

omega <- 1/(1 - c(G %*% lambda))/n_obs
# check: - sign difference in lambda from qin and lawless
omega - flexEL:::omega_hat(G)

sum(log(omega)) - logEL(G)

logEL2 <- function(G) {
  lambda <- lambdaNR(G)
  omega <- 1/(1 - c(G %*% lambda))/n_obs
  sum(log(omega))
}

matrix(numDeriv::grad(func = logEL2, x = G), n_obs, n_eq) -
  n_obs * outer(omega, lambda)

# ok more complex test

X0 <- matrix(rnorm(n_obs*n_eq), n_obs, n_eq)
y0 <- rnorm(n_obs)

Gfun <- function(theta) {
  mr_evalG(y0, X0, theta)
}

logEL3 <- function(theta) {
  G <- Gfun(theta)
  logEL2(G)
}

theta0 <- rnorm(n_eq)


G <- Gfun(theta0)
lambda <- lambdaNR(G)
omega <- 1/(1 - c(G %*% lambda))/n_obs
## dGdt <- array(jacobian(Gfun, x = theta0, X = X0), dim = c(n_obs, n_eq, n_eq))

dldt1 <- matrix(NA, n_obs, n_eq)
for(ii in 1:n_obs) {
  jac <- jacobian(function(theta) Gfun(theta)[ii,], x = theta0)
  dldt1[ii,] <- omega[ii] * c(jac %*% lambda)
}
dldt1 <- colSums(dldt1)

dldt2 <- grad(logEL3, x = theta0)

n_obs * dldt1 - dldt2
