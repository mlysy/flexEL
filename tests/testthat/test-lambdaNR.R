# ---- testing lambdaNR ----

source("test_utils.R")
source("el_rfuns.R")

# library(testthat) # not loaded automatically
context("lambdaNR")

ntest <- 50

# Non-censored case:
# checking R and C++ implementations are equal
test_that("no censoring: lambda_R == lambda_cpp", {
  for(ii in 1:ntest) {
    n <- sample(10:20,1)
    p <- sample(1:(n-2), 1)
    G <- matrix(rnorm(n*p),n,p) # randomly generated G
    max_iter <- sample(c(2, 10, 100), 1)
    rel_tol <- runif(1, 1e-6, 1e-5)
    lambda_cpp <- lambdaNR(G = G,
                           max_iter = max_iter, rel_tol = rel_tol, support = FALSE, 
                           verbose = FALSE)
    nrout <- lambdaNR_R(G = G,
                        max_iter = max_iter, rel_tol = rel_tol)
    lambda_R <- nrout$lambda
    expect_equal(lambda_R, lambda_cpp)
  }
})

# checking optimality of the solution from C++
test_that("no censoring: lambda_cpp is optimal", {
  for(ii in 1:ntest) {
    n <- sample(10:20,1)
    p <- sample(1:(n-2), 1)
    G <- matrix(rnorm(n*p),n,p) # randomly generated G
    max_iter <- sample(c(5, 10, 100), 1)
    rel_tol <- runif(1, 1e-6, 1e-5)
    lambda_cpp <- lambdaNR(G = G,
                           max_iter = max_iter, rel_tol = rel_tol, verbose = FALSE)
    # check optimality by optim_proj if converged
    if (!any(is.na(lambda_cpp))) {
      ocheck <- optim_proj(xsol = lambda_cpp,
                           fun = function(lambda) Qfun(lambda, G),
                           plot = FALSE)
      expect_lt(max_xdiff(ocheck),0.01)
    }
  }
})
