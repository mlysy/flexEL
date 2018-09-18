# ---- testing lambdaNR ----
# library(bayesEL) # always load the package (with library)
library(optimCheck)
source("el-utils.R")
source("el-rfuns.R")
source("el-model.R")

# library(testthat) # not loaded automatically
context("lambdaNR_cens")

ntest <- 50

# Censored case:
# checking R and C++ implementations are equal
test_that("under censoring: lambda.R == lambda.cpp", {
  for(ii in 1:ntest) {
    n <- sample(10:20,1)
    p <- sample(1:(n-2), 1)
    G <- matrix(rnorm(n*p),n,p) # randomly generated G
    max_iter <- sample(c(2, 10, 100), 1)
    rel_tol <- runif(1, 1e-6, 1e-5)
    weights <- abs(rnorm(n))
    weights <- weights / sum(weights) * n # sum(weights) == n
    lambda.cpp <- lambdaNR_cens(G = G, weights = weights, 
                                max_iter = max_iter, rel_tol = rel_tol, verbose = FALSE)
    nrout <- lambdaNRC_R(G = G, weights, max_iter = max_iter, rel_tol = rel_tol, verbose = FALSE)
    lambda.R <- nrout$lambda
    # check R and C++ are equal
    expect_equal(lambda.cpp, lambda.R)
  }
})

# checking optimality of the solution from C++
test_that("under censoring: lambda.cpp is optimal", {
  for(ii in 1:ntest) {
    n <- sample(10:20,1)
    p <- sample(1:(n-2), 1)
    G <- matrix(rnorm(n*p),n,p) # randomly generated G
    max_iter <- sample(c(5, 10, 100), 1)
    rel_tol <- runif(1, 1e-6, 1e-5)
    weights <- abs(rnorm(n))
    weights <- weights / sum(weights) * n # sum(weights) == n
    lambda.cpp <- lambdaNR_cens(G = G, weights = weights, 
                                max_iter = max_iter, rel_tol = rel_tol, verbose = FALSE)
    # check optimality by optim_proj if converged
    if (!any(is.na(lambda.cpp))) {
      ocheck <- optim_proj(xsol = lambda.cpp,
                           fun = function(lambda) QfunCens(lambda, G, weights),
                           plot = FALSE)
      expect_lt(max.xdiff(ocheck),0.01)
    }
  }
})
