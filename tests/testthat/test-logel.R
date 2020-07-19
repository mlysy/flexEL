# ---- testing logel implementation in R and C++ are equal ----

source("el_rfuns.R")

# library(testthat) # not loaded automatically
context("logEL")

ntest <- 50

# Non-censored case:
test_that("logel_cpp == logel_R no censoring, no support correction", {
  for(ii in 1:ntest) {
    n <- sample(10:20,1)
    p <- sample(1:(n-2), 1)
    max_iter <- sample(c(10, 100, 500), 1)
    rel_tol <- runif(1, 1e-6, 1e-5)
    G <- matrix(rnorm(n*p),n,p) # random G here
    support <- FALSE
    logopt_cpp <- flexEL::logEL(G = G, support = support,
                                max_iter = max_iter, rel_tol = rel_tol, abs_tol = 1e-3,
                                return_omega = FALSE, verbose = FALSE)
    omegahat_R <- omega_hat_R(G = G, adjust = support,
                              max_iter = max_iter, rel_tol = rel_tol, abs_tol = 1e-3, verbose = FALSE)
    logopt_R <- logEL_R(omegahat_R, adjust = support)
    if (!is.infinite(logopt_cpp) & !is.infinite(logopt_R)) {
      expect_equal(logopt_cpp, logopt_R, tolerance = 1e-4)
    }
  }
})

test_that("logel_cpp == logel_R no censoring, with support correction", {
  for(ii in 1:ntest) {
    n <- sample(10:20,1)
    p <- sample(1:(n-2), 1)
    max_iter <- sample(c(10, 100, 500), 1)
    rel_tol <- runif(1, 1e-6, 1e-5)
    abs_tol <- runif(1, 1e-5, 1e-3)
    G <- matrix(rnorm(n*p),n,p) # random G here
    support <- TRUE
    logopt_cpp <- flexEL::logEL(G = G, support = support, 
                                max_iter = max_iter, rel_tol = rel_tol, abs_tol = abs_tol, 
                                return_omega = FALSE, verbose = FALSE)
    omegahat_R <- omega_hat_R(G = adjG_R(G), adjust = support,
                              max_iter = max_iter, rel_tol = rel_tol, abs_tol = abs_tol, verbose = FALSE)
    logopt_R <- logEL_R(omegahat_R, adjust = support)
    if (!is.infinite(logopt_cpp) & !is.infinite(logopt_R)) {
      expect_equal(logopt_cpp, logopt_R, tolerance = 1e-4)
    }
  }
})

# Censored case:
test_that("logel_cpp == logel_R right-censored, no support correction", {
  for(ii in 1:ntest) {
    n <- sample(10:20,1)
    p <- sample(1:(n-2), 1)
    max_iter <- sample(c(10, 100, 500), 1)
    rel_tol <- runif(1, 1e-6, 1e-5)
    abs_tol <- runif(1, 1e-5, 1e-3)
    G <- matrix(rnorm(n*p), n, p)
    deltas <- rep(1,n)
    numcens <- sample(round(n/2),1)
    censinds <- sample(n,numcens)
    deltas[censinds] <- 0
    epsilons <- rnorm(n)
    support <- FALSE
    logopt_cpp <- flexEL::logEL(G = G, delta = deltas, eps = epsilons, support = support, 
                                max_iter = max_iter, rel_tol = rel_tol, abs_tol = abs_tol, verbose = FALSE)
    omegahat_R <- omega_hat_R(G = G, deltas = deltas, epsilons = epsilons, adjust = support,
                              max_iter = max_iter, rel_tol = rel_tol, abs_tol = abs_tol, verbose = FALSE)$omegas
    logopt_R <- logEL_R(omegas = omegahat_R, epsilons = epsilons, deltas = deltas, adjust = support)
    if (!is.infinite(logopt_cpp) & !is.infinite(logopt_R)) {
      expect_equal(logopt_cpp, logopt_R, tolerance = 1e-4)
    }
  }
})

# Censored case:
# failed <- 0
test_that("logel_cpp == logel_R right-censored, with support correction", {
  for(ii in 1:ntest) {
    n <- sample(10:20,1)
    p <- sample(1:(n-2), 1)
    max_iter <- sample(c(10, 100, 500), 1)
    rel_tol <- runif(1, 1e-6, 1e-5)
    abs_tol <- runif(1, 1e-5, 1e-3)
    G <- matrix(rnorm(n*p), n, p)
    deltas <- rep(1,n)
    numcens <- sample(round(n/2),1)
    censinds <- sample(n,numcens)
    deltas[censinds] <- 0
    epsilons <- rnorm(n)
    support <- TRUE
    logopt_cpp <- flexEL::logEL(G = G, delta = deltas, eps = epsilons, support = support, 
                                max_iter = max_iter, rel_tol = rel_tol, abs_tol = abs_tol, verbose = FALSE)
    omegahat_R <- omega_hat_R(G = adjG_R(G), deltas = deltas, epsilons = epsilons, adjust = support,
                              max_iter = max_iter, rel_tol = rel_tol, abs_tol = abs_tol, verbose = FALSE)$omegas
    logopt_R <- logEL_R(omegas = omegahat_R, epsilons = epsilons, deltas = deltas, adjust = support)
    if (!is.infinite(logopt_cpp) & !is.infinite(logopt_R)) {
      expect_equal(logopt_cpp, logopt_R, tolerance = 1e-4)
    }
    # else {
    #   failed <<- failed + 1
    #   message("R version did not converge")
    # }
  }
})

