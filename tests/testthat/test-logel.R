# ---- testing logel implementation in R and C++ are equal ----
## library(bayesEL) # always load the package (with library)
# library(optimCheck)
source("el-utils.R")
## source("~/bayesEL/tests/testthat/el-utils.R")

# library(testthat) # not loaded automatically
context("logEL")

ntest <- 50

# Non-censored case:
test_that("logel.R == logel.cpp", {
  for(ii in 1:ntest) {
    n <- sample(10:20,1)
    p <- sample(1:(n-5), 1)
    max_iter <- sample(c(2, 10, 100), 1)
    rel_tol <- runif(1, 1e-6, 1e-5)
    G <- matrix(rnorm(n*p),n,p)
    logopt.cpp <- logEL(G,max_iter=max_iter,rel_tol=rel_tol)
    logopt.R <- logEL_R(G,max_iter=max_iter,rel_tol=rel_tol)
    expect_equal(logopt.cpp,logopt.R)
  }
})

# Censored case:
test_that("logelC.R == logelC.cpp", {
  for(ii in 1:ntest) {
    n <- sample(10:20,1)
    p <- sample(1:(n-5), 1)
    max_iter <- sample(c(2, 10, 100), 1)
    rel_tol <- runif(1, 1e-6, 1e-5)
    G <- matrix(rnorm(n*p),n,p)
    deltas <- rep(1,n)
    numcens <- sample(round(n/2),1)
    censinds <- sample(n,numcens)
    deltas[censinds] <- 0
    epsilons <- rnorm(n)
    logopt.cpp <- logEL(G,deltas,epsilons,max_iter=max_iter,rel_tol=rel_tol)
    # logopt.cpp
    logopt.R <- logEL_R(G,deltas,epsilons,max_iter=max_iter,rel_tol=rel_tol)
    # logopt.R
    if ((logopt.R == -Inf) && (logopt.cpp != -Inf)) {
      message("R version did not converge but C++ did.")
    }
    else {
      expect_equal(logopt.cpp,logopt.R)
    }
    # location model data generation
    # X <- matrix(rnorm(n*p),n,p)
    # beta <- 2+ rnorm(p)
    # y <- c(X %*% beta + rnorm(n))
    # G <- mr.evalG(y,X,beta)
    # Old code: find the maximized loglikelihood
    # logEL(G, max_iter = max_iter, rel_tol = rel_tol, verbose = TRUE)
  }
})
