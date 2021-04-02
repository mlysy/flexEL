library(testthat)
context("CensEL")

source("el_rfuns.R")

ntest <- 5

# ---- eval_weights ----

test_that("eval_weights with default settings", {
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
    omega <- runif(n, 0, 1)
    omega <- omega/sum(omega)
    weights_cpp <- cel$eval_weights(delta, epsilon, omega)
    weights_R <- evalWeights_R(delta, omega, epsilon)
    # flexEL::.EvalWeights(omega, delta, epsilon, FALSE)
    expect_equal(weights_cpp, weights_R)
  }
})

# ---- omega_hat -----

# nconv <- 0
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
    omega_cpp
    omega_R_lst <- omega_hat_R(G, deltas = delta, epsilons = epsilon)
    omega_R_lst$omegas
    range(omega_cpp-omega_R_lst$omegas)
    # omega_R_lst$conv
    if (!all(is.na(omega_cpp)) & omega_R_lst$conv) {
      # nconv <<- nconv + 1
      expect_equal(omega_cpp, omega_R_lst$omegas, tolerance = 1e-2)
    }
  }
})

# TODO: discrepancy too obvious .. need to check algorithm
# nconv <- 0
test_that("omega_hat with given convergence settings", {
  for (ii in 1:ntest) {
    n <- sample(10:20,1)
    p <- sample(1:(n-2), 1)
    max_iter <- sample(c(100, 200, 500), 1)
    rel_tol <- runif(1, 1e-7, 1e-3)
    abs_tol <- runif(1, 1e-4, 1e-2)
    cel <- CensEL$new(n, p)
    cel$max_iter_nr <- max_iter
    cel$rel_tol <- rel_tol
    cel$max_iter_em <- max_iter
    cel$abs_tol <- abs_tol
    G <- matrix(rnorm(n*p),n,p)
    delta <- rep(1,n)
    numcens <- sample(round(n/2),1)
    censinds <- sample(n,numcens)
    delta[censinds] <- 0
    epsilon <- rnorm(n)
    omega_cpp <- cel$omega_hat(G, delta = delta, epsilon = epsilon)
    # omega_cpp
    omega_R_lst <- omega_hat_R(G, deltas = delta, epsilons = epsilon, 
                               max_iter = max_iter, rel_tol = rel_tol, abs_tol = abs_tol)
    # omega_R_lst$omegas
    range(omega_cpp-omega_R_lst$omegas)
    # omega_R_lst$conv
    if (!all(is.na(omega_cpp)) & omega_R_lst$conv) {
      # nconv <<- nconv + 1
      expect_equal(omega_cpp, omega_R_lst$omegas, tolerance = 0.05)
    }
  }
})

# TODO: discrepancy too obvious .. need to check algorithm
nconv <- 0
test_that("omega_hat with given convergence settings and support correction", {
  for (ii in 1:ntest) {
    n <- sample(10:20,1)
    p <- sample(1:(n-2), 1)
    max_iter <- sample(c(100, 200, 500), 1)
    rel_tol <- runif(1, 1e-7, 1e-4)
    abs_tol <- runif(1, 1e-4, 1e-2)
    adj_a <- runif(1, 1, 5)
    cel <- CensEL$new(n, p)
    cel$max_iter_nr <- max_iter
    cel$rel_tol <- rel_tol
    cel$max_iter_em <- max_iter
    cel$abs_tol <- abs_tol
    cel$supp_adj <- TRUE
    cel$supp_adj_a <- adj_a
    G <- matrix(rnorm(n*p),n,p)
    delta <- rep(1,n)
    numcens <- sample(round(n/2),1)
    censinds <- sample(n,numcens)
    delta[censinds] <- 0
    epsilon <- rnorm(n)
    omega_cpp <- cel$omega_hat(G, delta = delta, epsilon = epsilon)
    omega_cpp
    omega_R_lst <- omega_hat_EM_R(adjG_R(G, adj_a), deltas = delta, epsilons = epsilon, adjust = TRUE,
                               max_iter = max_iter, rel_tol = rel_tol, abs_tol = abs_tol)
    omega_R_lst$omegas
    omega_R_lst$conv
    if (!all(is.na(omega_cpp)) & omega_R_lst$conv) {
      nconv <<- nconv + 1
      expect_equal(omega_cpp, omega_R_lst$omegas, tolerance = abs_tol)
    }
  }
})

