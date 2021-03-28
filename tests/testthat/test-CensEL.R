library(testthat)
context("CensEL")

source("el_rfuns.R")

ntest <- 5

# ---- eval_weights ----

n <- sample(10:20,1)
p <- sample(1:(n-2), 1)
cel <- CensEL$new(n, p)
G <- matrix(rnorm(n*p),n,p)
delta <- rep(1,n)
numcens <- sample(round(n/2),1)
censinds <- sample(n,numcens)
delta[censinds] <- 0
epsilon <- rnorm(n)
omega <- runif(n, 0, 1)
omega <- omega/sum(omega)
weights_cpp <- cel$eval_weights(delta, epsilon, omega)

# ---- omega_hat -----

nconv <- 0
test_that("omega_hat with default settings", {
  for (ii in 1:ntest) {
    n <- sample(10:20,1)
    p <- sample(1:(n-2), 1)
    cel <- CensEL$new(n, p)
    G <- matrix(rnorm(n*p),n,p)
    delta <- rep(1,n)
    numcens <- sample(round(n/2),1)
    censinds <- sample(n,numcens)
    delta[censinds] <- 0
    epsilon <- rnorm(n)
    omega_cpp <- cel$omega_hat(G, delta = delta, epsilon = epsilon)
    omega_R <- omega_hat_R(G)
    if (!any(is.na(omega_cpp)) & !any(is.na(omega_R))) {
      # nconv <<- nconv + 1
      expect_equal(omega_cpp, omega_R)
    }
    else {
      expect_equal(any(is.na(omega_cpp)), any(is.na(omega_R)))
    }
  }
})
