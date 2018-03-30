#---- testing omega.hat  ----
library(bayesEL) # always load the package (with library)
library(optimCheck)
# source("el-utils.R")
source("~/bayesEL/tests/testthat/el-utils.R")

# library(testthat) # not loaded automatically
context("omega.hat")

ntest <- 50

# Non-censored case: 
# checking R and C++ implementations are equal (not optimality) 
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
# TODO: check the correctness of the R version by the old 5 variables way (?)
# and then check the c++ version against it 
# test_that("omegahatC.R == omegahatC.cpp", {
#     for(ii in 1:ntest) {
#         # TODO: may should choose a random percentage of censored observations
#         # c <- abs(rnorm(n, mean=y, sd=1)) # independent random censoring values
#         # deltas <- y <= c # delta == 1 if observed 0 if censored
#         # y[!deltas] <- c[!deltas] # observed lifetime after censoring
#         # 1-sum(deltas)/n
#         # debug:
#         n <- 15
#         p <- 2
#         X <- replicate(p, rnorm(n))
#         X[1,] <- rep(1,p)
#         beta0 <- rnorm(p)
#         y <- c(X %*% beta0) + rnorm(n) # with N(0,1) error term
#         max_iter <- 100
#         rel_tol <- 1e-5
#         epsilons <- c(y - X %*% beta0)
#         deltas <- rep(1,n)
#         # censind <- sample(n,2)
#         censind <- order(epsilons)[c(4,7)]
#         deltas[censind] <- 0
#         y[censind] <- rnorm(1,mean=y[censind]) # makes only 1 censored
#         G.cpp <- mr.evalG(y, X, beta0)
#         omegahat.cpp <- omega.hat(G.cpp, deltas, epsilons, max_iter = max_iter, rel_tol = rel_tol, verbose = FALSE)
#         omegahat.R <- omega.hat_R(G.R, deltas, epsilons, max_iter = max_iter, rel_tol = rel_tol, verbose = FALSE)
#         omegahat.cpp
#         omegahat.R
#         omegahat.cpp - omegahat.R
#     }
# })
