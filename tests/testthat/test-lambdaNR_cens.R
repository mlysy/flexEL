# ---- testing lambdaNR ----

source("test_utils.R")
source("el_rfuns.R")

# library(testthat) # not loaded automatically
context("lambdaNR_cens")

ntest <- 50

# Censored case:
# checking R and C++ implementations are equal
test_that("under censoring: lambda_R == lambda_cpp", {
  for(ii in 1:ntest) {
    n <- sample(10:20,1)
    p <- sample(1:(n-2), 1)
    G <- matrix(rnorm(n*p),n,p) # randomly generated G
    max_iter <- sample(c(2, 10, 100), 1)
    rel_tol <- runif(1, 1e-6, 1e-5)
    weights <- abs(rnorm(n))
    weights <- weights / sum(weights) * n # sum(weights) == n
    lambda_cpp <- lambdaNR_cens(G = G, weights = weights, 
                                max_iter = max_iter, rel_tol = rel_tol, verbose = FALSE)
    nrout <- lambdaNRC_R(G = G, weights, max_iter = max_iter, rel_tol = rel_tol, verbose = FALSE)
    lambda_R <- nrout$lambda
    # check R and C++ are equal
    expect_equal(lambda_cpp, lambda_R)
  }
})

# checking optimality of the solution from C++
test_that("under censoring: lambda_cpp is optimal", {
  for(ii in 1:ntest) {
    n <- sample(10:20,1)
    p <- sample(1:(n-2), 1)
    G <- matrix(rnorm(n*p),n,p) # randomly generated G
    max_iter <- sample(c(5, 10, 100), 1)
    rel_tol <- runif(1, 1e-6, 1e-5)
    weights <- abs(rnorm(n))
    weights <- weights / sum(weights) * n # sum(weights) == n
    lambda_cpp <- lambdaNR_cens(G = G, weights = weights, 
                                max_iter = max_iter, rel_tol = rel_tol, verbose = FALSE)
    # check optimality by optim_proj if converged
    if (!any(is.na(lambda_cpp))) {
      ocheck <- optim_proj(xsol = lambda_cpp,
                           fun = function(lambda) QfunCens(lambda, G, weights),
                           plot = FALSE)
      expect_lt(max_xdiff(ocheck),0.01)
    }
  }
})
