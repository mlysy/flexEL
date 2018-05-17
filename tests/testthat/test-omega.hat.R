# ---- testing R and C++ implementations of omega.hat are equal and optimality ----
# library(bayesEL) # always load the package (with library)
library(optimCheck)
source("el-utils.R")
source("el-rfuns.R")
source("el-model.R")

# library(testthat) # not loaded automatically
context("omega.hat")

ntest <- 50

# Non-censored case:
# checking R and C++ implementations are equal
test_that("no censoring: omegahat.R == omegahat.cpp", {
  for(ii in 1:ntest) {
    n <- sample(10:20,1)
    p <- sample(1:(n-2), 1)
    # max_iter <- sample(c(2, 10, 100), 1)
    max_iter <- sample(c(10, 100, 500), 1)
    rel_tol <- runif(1, 1e-6, 1e-5)
    G <- matrix(rnorm(n*p),n,p) # random G here
    omegahat.cpp <- omega.hat(G = G, max_iter = max_iter, rel_tol = rel_tol, verbose = FALSE)
    omegahat.R <- omega.hat_R(G = G, max_iter = max_iter, rel_tol = rel_tol, verbose = FALSE)
    expect_equal(omegahat.cpp, omegahat.R)
  }
})

# checking optimality of the solution from C++
test_that("no censoring: omegahat.cpp is optimal", {
  for(ii in 1:ntest) {
    n <- sample(10:20,1)
    p <- sample(1:(n-2), 1)
    # max_iter <- sample(c(5, 10, 100), 1)
    max_iter <- sample(c(10, 100, 500), 1)
    rel_tol <- runif(1, 1e-6, 1e-5)
    G <- matrix(rnorm(n*p),n,p) # random G here
    omegahat.cpp <- omega.hat(G = G, max_iter = max_iter, rel_tol = rel_tol, verbose = FALSE)
    # check optimality by optim_proj if converged
    if (!any(is.nan(omegahat.cpp))) {
      # optimal of omega.check should be at xsol = 1s if omegahat is optimal
      ocheck <- optim_proj(xsol = rep(1,n-p),
                           xrng = 0.05,
                           fun = function(x) {omega.check(x, omegahat.cpp, G)},
                           plot = FALSE)
      expect_lt(max.xdiff(ocheck),0.01)
    }
  }
})

# Censored case:
test_that("under censoring: omegahatC.R == omegahatC.cpp", {
  for(ii in 1:ntest) {
    n <- sample(10:20,1)
    p <- sample(1:(n-2), 1)
    # max_iter <- sample(c(2, 10, 100), 1)
    max_iter <- sample(c(10, 100, 500), 1)
    rel_tol <- runif(1, 1e-6, 1e-5)
    G <- matrix(rnorm(n*p), n, p)
    deltas <- rep(1,n)
    numcens <- sample(round(n/2),1)
    censinds <- sample(n,numcens)
    deltas[censinds] <- 0
    epsilons <- rnorm(n)
    omegahat.cpp <- omega.hat(G, deltas, epsilons, max_iter = max_iter, rel_tol = rel_tol, verbose = FALSE)
    omegahat.R <- omega.hat_R(G, deltas, epsilons, max_iter = max_iter, rel_tol = rel_tol, verbose = FALSE)
    if (!any(is.nan(omegahat.cpp)) && any(is.nan(omegahat.R))) {
      message("R version did not converge but C++ does.")
    }
    else {
      expect_equal(omegahat.cpp, omegahat.R)
    }
  }
})

# checking optimality of the solution from C++
test_that("under censoring: omegahat.cpp is optimal", {
  for(ii in 1:ntest) {
    n <- sample(10:20,1)
    p <- sample(1:(n-2), 1)
    # max_iter <- sample(c(2, 10, 100), 1)
    max_iter <- sample(c(10, 100, 500), 1)
    rel_tol <- runif(1, 1e-6, 1e-5)
    G <- matrix(rnorm(n*p), n, p)
    deltas <- rep(1,n)
    numcens <- sample(round(n/2),1)
    censinds <- sample(n,numcens)
    deltas[censinds] <- 0
    epsilons <- rnorm(n)
    omegahat.cpp <- omega.hat(G, deltas, epsilons, max_iter = max_iter, rel_tol = rel_tol)
  }
  if (!any(is.nan(omegahat.cpp))) {
    ocheck <- optim_proj(xsol = rep(1,n-p),
                         xrng = 0.05,
                         npts = 101, # 101 would contain x exactly 1, o.w. sometimes does not work
                         fun = function(x) {omega.check(x, omegahat.cpp, G, deltas, epsilons)},
                         plot = FALSE)
    # print(ocheck)
    expect_lt(max.xdiff(ocheck), 0.01)
  }
})


## omega.hat <- function(G, deltas, lambda) {
## }

## # method 1:
## omega <- omega.hat(G, deltas, ...)
## # method 2:
## lambda <- lambdaNR(G, deltas, ...)
## omega <- omega.hat(G, lambda)
