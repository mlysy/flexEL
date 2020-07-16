# ---- testing logel implementation in R and C++ are equal ----

source("el_rfuns.R")

# library(testthat) # not loaded automatically
context("logEL")

ntest <- 50

# Non-censored case:
test_that("logel_R == logel_cpp no censoring, no support correction", {
  for(ii in 1:ntest) {
    n <- sample(10:20,1)
    p <- sample(1:(n-2), 1)
    max_iter <- sample(c(10, 100, 500), 1)
    rel_tol <- runif(1, 1e-6, 1e-5)
    G <- matrix(rnorm(n*p),n,p) # random G here
    support <- FALSE
    logopt_cpp <- logEL(G = G, support = support, 
                        max_iter = max_iter, rel_tol = rel_tol, abs_tol = 1e-3, 
                        return_omega = FALSE, verbose = FALSE)
    omegahat_R <- omega_hat_R(G = G, adjust = support,
                              max_iter = max_iter, rel_tol = rel_tol, abs_tol = 1e-3, verbose = FALSE)
    logopt_R <- logEL_R(omegahat_R, adjust = support)
    expect_equal(logopt_cpp, logopt_R)
  }
})

test_that("logel_R == logel_cpp no censoring, with support correction", {
  for(ii in 1:ntest) {
    n <- sample(10:20,1)
    p <- sample(1:(n-2), 1)
    max_iter <- sample(c(10, 100, 500), 1)
    rel_tol <- runif(1, 1e-6, 1e-5)
    G <- matrix(rnorm(n*p),n,p) # random G here
    support <- TRUE
    logopt_cpp <- logEL(G = G, support = support, 
                        max_iter = max_iter, rel_tol = rel_tol, abs_tol = 1e-3, 
                        return_omega = FALSE, verbose = FALSE)
    omegahat_R <- omega_hat_R(G = adjG_R(G), adjust = support,
                              max_iter = max_iter, rel_tol = rel_tol, abs_tol = 1e-3, verbose = FALSE)
    logopt_R <- logEL_R(omegahat_R, adjust = support)
    expect_equal(logopt_cpp, logopt_R)
  }
})

# Censored case:
test_that("logel_R == logel_cpp right-censored, no support correction", {
  for(ii in 1:ntest) {
    n <- sample(10:20,1)
    omegas <- abs(rnorm(n))
    k <- sample(round(n/2),1)
    idx <- sample(n,k)
    omegas[idx] <- omegas[idx]/(10^sample(5:25,k))
    omegas <- omegas/sum(omegas)
    deltas <- rep(1,n)
    deltas[idx] <- 0
    epsilons <- rnorm(n)
    logopt_cpp <- logEL(omegas,epsilons,deltas)
    logopt_R <- logEL_R(omegas,epsilons,deltas)
    expect_equal(logopt_cpp,logopt_R)
  }
})
