# ---- testing R and C++ implementations of logEL_smooth are equal ----

source("el_rfuns.R")

# library(testthat) # not loaded automatically
context("logEL_smooth")

ntest <- 50

test_that("logEL_smooth_R == logEL_smooth_cpp, no support correction", {
  
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
    sp <- sample(1:100, 1)
    logopt_cpp <- flexEL::logEL(G = G, delta = deltas, eps = epsilons, support = support, sp = sp,
                                max_iter = max_iter, rel_tol = rel_tol, abs_tol = abs_tol, verbose = FALSE)
    omegahat_R <- omega_hat_R(G = G, deltas = deltas, epsilons = epsilons, adjust = support,
                              max_iter = max_iter, rel_tol = rel_tol, abs_tol = abs_tol, verbose = FALSE)$omegas
    logopt_R <- logEL_smooth_R(omegas = omegahat_R, epsilons = epsilons, deltas = deltas, s = sp, adjust = support)
    if (!is.infinite(logopt_cpp) & !is.infinite(logopt_R)) {
      expect_equal(logopt_cpp, logopt_R, tolerance = 1e-4)
    }
    # else {
    #   failed <<- failed + 1
    #   message("R version did not converge")
    # }
  }
})


test_that("logEL_smooth_R == logEL_smooth_cpp, with support correction", {
  
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
    sp <- sample(1:100, 1)
    logopt_cpp <- flexEL::logEL(G = G, delta = deltas, eps = epsilons, support = support, sp = sp,
                                max_iter = max_iter, rel_tol = rel_tol, abs_tol = abs_tol, verbose = FALSE)
    omegahat_R <- omega_hat_R(G = adjG_R(G), deltas = deltas, epsilons = epsilons, adjust = support,
                              max_iter = max_iter, rel_tol = rel_tol, abs_tol = abs_tol, verbose = FALSE)$omegas
    logopt_R <- logEL_smooth_R(omegas = omegahat_R, epsilons = epsilons, deltas = deltas, s = sp, adjust = support)
    if (!any(is.nan(omegahat_R))) {
      expect_equal(logopt_cpp, logopt_R, tolerance = 1e-4)
    }
    # else {
    #   failed <<- failed + 1
    #   message("R version did not converge")
    # }
  }
})
