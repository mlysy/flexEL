
library(testthat)
source("el_rfuns.R")

ntest <- 10

# nconv <- 0
test_that("lambda_nr with default options", {
  for (ii in 1:ntest) {
    n <- sample(10:20,1)
    p <- sample(1:(n-2), 1)
    gel <- GenEL$new(n, p)
    G <- matrix(rnorm(n*p),n,p)
    lambda_cpp <- gel$lambda_nr(G)
    lambda_R_lst <- lambdaNR_R(G)
    if (!lambda_R_lst$convergence) {
      expect_equal(all(is.na(lambda_cpp)), all(is.na(lambda_R_lst$lambda)))
    }
    else {
      # nconv <<- nconv + 1
      expect_equal(lambda_cpp, lambda_R_lst$lambda)
    }
  }
})

# nconv <- 0
test_that("lambda_nr with given convergence settings", {
  for (ii in 1:ntest) {
    n <- sample(10:20,1)
    p <- sample(1:(n-2), 1)
    max_iter <- sample(c(100, 500), 1)
    rel_tol <- runif(1, 1e-6, 1e-2)
    gel <- GenEL$new(n, p)
    gel$max_iter <- max_iter
    gel$rel_tol <- rel_tol
    G <- matrix(rnorm(n*p),n,p)
    lambda_cpp <- gel$lambda_nr(G, verbose = TRUE)
    lambda_R_lst <- lambdaNR_R(G, max_iter = max_iter, rel_tol = rel_tol)
    if (!lambda_R_lst$convergence) {
      expect_equal(all(is.na(lambda_cpp)), all(is.na(lambda_R_lst$lambda)))
    }
    else {
      # nconv <<- nconv + 1
      expect_equal(lambda_cpp, lambda_R_lst$lambda)
    }
  }
})

# nconv <- 0
test_that("omega_hat with default settings", {
  for (ii in 1:ntest) {
    n <- sample(10:20,1)
    p <- sample(1:(n-2), 1)
    gel <- GenEL$new(n, p)
    G <- matrix(rnorm(n*p),n,p)
    omega_cpp <- gel$omega_hat(G)
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

