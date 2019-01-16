# ---- testing R and C++ implementations of omega.hat are equal and optimality ----
# library(bayesEL) # always load the package (with library)
library(optimCheck)
source("el-utils.R")
source("el-rfuns.R")
source("el-model.R")
source("../dontrun/smoothEL.R")

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

# checking optimality of the solution from C++ (with adjusted G matrix) 
# (TODO: make adjG internal..)
test_that("no censoring: omegahat.cpp is optimal", {
  for(ii in 1:ntest) {
    n <- sample(10:20,1)
    p <- sample(1:(n-2), 1)
    # max_iter <- sample(c(5, 10, 100), 1)
    max_iter <- sample(c(10, 100, 500), 1)
    rel_tol <- runif(1, 1e-6, 1e-5)
    G <- matrix(rnorm(n*p),n,p) # random G here
    G <- adjG_R(G)
    omegahat.cpp <- omega.hat(G = G, max_iter = max_iter, rel_tol = rel_tol, verbose = FALSE)
    ocheck <- optim_proj(xsol = rep(1,n-p+1),
                         xrng = 0.05,
                         fun = function(x) {omega.check(x, omegahat.cpp, G)},
                         plot = FALSE)
    expect_lt(max.xdiff(ocheck),0.01)
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
    abs_tol <- runif(1, 1e-5, 1e-3)
    G <- matrix(rnorm(n*p), n, p)
    deltas <- rep(1,n)
    numcens <- sample(round(n/2),1)
    censinds <- sample(n,numcens)
    deltas[censinds] <- 0
    epsilons <- rnorm(n)
    omegahat.cpp <- omega.hat(G, deltas, epsilons, max_iter = max_iter, 
                              rel_tol = rel_tol, abs_tol = abs_tol, verbose = FALSE)
    omegahat.R <- omega.hat_R(G, deltas, epsilons, max_iter = max_iter, 
                              rel_tol = rel_tol, abs_tol = abs_tol, verbose = FALSE)$omegas
    if (!any(is.nan(omegahat.cpp)) && any(is.nan(omegahat.R))) {
      message("R version did not converge but C++ does.")
    }
    else {
      expect_equal(omegahat.cpp, omegahat.R)
    }
  }
})

# save.ome <- list()
# save.del <- list()
# save.G <- list()
# save.eps <- list()
# save.ite <- list()
# save.tol <- list()
# count <- 0

# checking optimality of the solution from C++
test_that("under censoring: omegahat.cpp is optimal", {
  for(ii in 1:ntest) {
    n <- sample(10:30,1)
    p <- sample(1:(n-2), 1)
    # max_iter <- sample(c(2, 10, 100), 1)
    max_iter <- sample(c(10, 100, 500), 1)
    # max_iter <- 200
    rel_tol <- runif(1, 1e-8, 1e-6)
    abs_tol <- runif(1, 1e-5, 1e-3)
    # rel_tol <- 1e-7
    G <- matrix(rnorm(n*p), n, p)
    deltas <- rep(1,n)
    numcens <- sample(round(n/2),1)
    censinds <- sample(n,numcens)
    deltas[censinds] <- 0
    epsilons <- rnorm(n)
    omegahat.cpp <- omega.hat(G, deltas, epsilons, max_iter = max_iter, 
                              rel_tol = rel_tol, abs_tol = abs_tol, verbose = FALSE)
    if (!any(is.nan(omegahat.cpp))) {
      # ocheck <- optim_proj(xsol = rep(1,n-p),
      #                      xrng = 0.01,
      #                      npts = 501, # 101 would contain x exactly 1, o.w. sometimes does not work
      #                      fun = function(x) {omega.check(x, omegahat.cpp, G, deltas, epsilons)},
      #                      plot = FALSE)
      # count <- count + 1
      # save.ome[[count]] <- omegahat.cpp
      # save.del[[count]] <- deltas
      # save.G[[count]] <- G
      # save.eps[[count]] <- epsilons
      # save.ite[[count]] <- max_iter
      # save.tol[[count]] <- rel_tol
      idx0 <- (abs(omegahat.cpp) < 1e-5 & !deltas)
      if (n-p-sum(idx0) > 0) {
        ocheck <- optim_proj(xsol = rep(1,n-p-sum(idx0)),
                             xrng = 0.01,
                             npts = 201, 
                             fun = function(x) {omega.pcheck(x, omegahat.cpp, G, deltas, epsilons, idx0, rel_tol)},
                             plot = FALSE)
        # print(ocheck)
        expect_lt(max.xdiff(ocheck), 0.01)
      }
    }
  }
})

# smoothed censored EL
# checking optimality of the solution (TODO: the following is in R)
# n <- 25
# p <- 5
# max_iter <- 200
# rel_tol <- 1e-5
# G <- matrix(rnorm(n*p), n, p)
# deltas <- rep(1,n)
# numcens <- sample(round(n/3),1)
# censinds <- sample(n,numcens)
# deltas[censinds] <- 0
# 1-sum(deltas)/n
# epsilons <- rnorm(n)
# oout <- omega.hat.EM.smooth_R(G, deltas, epsilons, max_iter = max_iter, rel_tol = rel_tol, verbose = FALSE)
# oout$conv
# omegahat <- oout$omegas
# ocheck <- optim_proj(xsol = rep(1,n-p),
#                      xrng = 0.0001,
#                      npts = 201,
#                      fun = function(x) {omega.smooth.check(x, omegahat, G, deltas, epsilons)},
#                      plot = FALSE)
# expect_lt(max.xdiff(ocheck),0.01)

# Censored case + smooth:
test_that("under censoring: omegahatCS.R == omegahatCS.cpp", {
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
    s <- sample(1:100,1)
    support <- FALSE
    omegahat.cpp <- omega.hat.EM.smooth(G, deltas, epsilons, s, 
                                        max_iter = max_iter, rel_tol = rel_tol, abs_tol = rel_tol, 
                                        support = support, verbose = FALSE)
    omegahat.R <- omega.hat.EM.smooth_R(G, deltas, epsilons, s, max_iter = max_iter, rel_tol = rel_tol, abs_tol = rel_tol, verbose = FALSE)$omegas
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
    n <- sample(10:30,1)
    p <- sample(1:(n-2), 1)
    # max_iter <- sample(c(2, 10, 100), 1)
    max_iter <- sample(c(10, 100, 500), 1)
    # max_iter <- 200
    rel_tol <- runif(1, 1e-8, 1e-6)
    # rel_tol <- 1e-7
    G <- matrix(rnorm(n*p), n, p)
    deltas <- rep(1,n)
    numcens <- sample(round(n/4),1)
    censinds <- sample(n,numcens)
    deltas[censinds] <- 0
    1-sum(deltas)/n
    epsilons <- rnorm(n)
    s <- sample(1:100,1)
    support <- FALSE
    omegahat.cpp <- omega.hat.EM.smooth(G, deltas, epsilons, s, max_iter = max_iter, 
                                          rel_tol = 1e-5, abs_tol = 1e-5, support = support, verbose = FALSE)
    # omegahat.cpp <- omega.hat.EM.smooth_R(G, deltas, epsilons, s, max_iter = max_iter,
    #                                       rel_tol = 1e-3, verbose = FALSE)$omegas
    if (!any(is.nan(omegahat.cpp))) {
      idx0 <- (abs(omegahat.cpp) < 1e-5 & !deltas)
      if (n-p-sum(idx0) > 0) {
        ocheck <- optim_proj(xsol = rep(1,n-p-sum(idx0)),
                             xrng = 0.01,
                             npts = 201,
                             fun = function(x) {omega.smooth.pcheck(x, omegahat.cpp, G, deltas, epsilons, idx0, s)},
                             plot = FALSE)
        # ocheck <- optim_proj(xsol = rep(1,n-p),
        #                      xrng = 0.001,
        #                      npts = 201, 
        #                      fun = function(x) {omega.smooth.check(x, omegahat.cpp, G, deltas, epsilons, s)},
        #                      plot = FALSE)
        # print(ocheck)
        expect_lt(max.xdiff(ocheck), 0.01)
      }
    }
  }
})

# --------------------------------------------------------------------------- 

## omega.hat <- function(G, deltas, lambda) {
## }

## # method 1:
## omega <- omega.hat(G, deltas, ...)
## # method 2:
## lambda <- lambdaNR(G, deltas, ...)
## omega <- omega.hat(G, lambda)
