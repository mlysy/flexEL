#---- testing omega.hat  ----
library(bayesEL) # always load the package (with library)
library(optimCheck)
# source("el-utils.R")
source("~/bayesEL/tests/testthat/el-utils.R")

# library(testthat) # not loaded automatically
context("omega.hat")

ntest <- 50

# Non-censored case: 
# checking R and C++ implementations are equal and optimality of omega.hat
test_that("omegahat.R == omegahat.cpp", {
  for(ii in 1:ntest) {
    n <- sample(10:20,1)
    p <- sample(1:(n-2), 1)
    max_iter <- sample(c(2, 10, 100), 1)
    rel_tol <- runif(1, 1e-6, 1e-5)
    G <- matrix(rnorm(n*p),n,p) # random G here
    omegahat.cpp <- omega.hat(G = G, max_iter = max_iter, rel_tol = rel_tol, verbose = FALSE)
    omegahat.R <- omega.hat_R(G = G, max_iter = max_iter, rel_tol = rel_tol, verbose = FALSE)
    # omegahat.cpp
    expect_equal(omegahat.cpp, omegahat.R)
    # optimal of omega.check should be at xsol = 1s if omegahat is optimal
    if (sum(omegahat.cpp)!=0) {
      ocheck <- optim_proj(xsol = rep(1,n-p), 
                           xrng = 0.05,
                           fun = function(x) {omega.check(x, omegahat.cpp, G)},
                           plot = FALSE) 
      expect_lt(max.xdiff(ocheck),0.01)
    }
  }
})

# Censored case: 
test_that("omegahatC.R == omegahatC.cpp", {
  for(ii in 1:ntest) {
    # message(ii)
    n <- sample(10:20,1)
    p <- sample(1:(n-2), 1)
    max_iter <- sample(c(2, 10, 100), 1)
    rel_tol <- runif(1, 1e-6, 1e-5)
    G <- matrix(rnorm(n*p), n, p)
    deltas <- rep(1,n)
    numcens <- sample(round(n/2),1)
    censinds <- sample(n,numcens)
    deltas[censinds] <- 0
    # deltas
    epsilons <- rnorm(n)
    omegahat.cpp <- omega.hat(G, deltas, epsilons, max_iter = max_iter, rel_tol = rel_tol, verbose = FALSE)
    # omegahat.cpp
    omegahat.R <- omega.hat_R(G, deltas, epsilons, max_iter = max_iter, rel_tol = rel_tol, verbose = FALSE)
    # omegahat.R
    if (sum(omegahat.cpp)!=0 && sum(omegahat.R)==0) {
      message("R version did not converge but C++ does, still checks optimality!")
    }
    else {
      expect_equal(omegahat.cpp, omegahat.R)
    }
    if (sum(omegahat.cpp)!=0) {
      ocheck <- optim_proj(xsol = rep(1,n-p), 
                           xrng = 0.05,
                           npts = 101, # 101 would have x exactly 1, o.w. sometimes does not work
                           fun = function(x) {omega.check(x, omegahat.cpp, G, deltas, epsilons)},
                           plot = FALSE) 
      # print(ocheck)
      expect_lt(max.xdiff(ocheck), 0.01)
    }
  }
})
