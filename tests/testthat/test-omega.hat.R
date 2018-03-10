#---- testing omega.hat  ----
library(bayesEL) # always load the package (with library)
# source("el-utils.R")
source("~/bayesEL/tests/testthat/el-utils.R")
source("~/bayesEL/tests/testthat/mle-check.R")

# library(testthat) # not loaded automatically
context("omega.hat")

ntest <- 50

test_that("lambda.R == lambda.cpp", {
    for(ii in 1:ntest) {
        # n <- sample(10:20,1)
        # p <- sample(1:(n-2), 1)
        # G <- matrix(rnorm(n*p),n,p) # random G seems not easy to work 
        system.time({
        # n <- sample(10:20,1)
        # p <- sample(1:(n-2), 1)
        n <- 15
        p <- 2
        X <- replicate(p, rnorm(n))
        X[1,] <- rep(1,p)
        beta0 <- rnorm(p)
        y <- c(X %*% beta0) + rnorm(n) # with N(0,1) error term
        lambda0 <- rnorm(p) # m == p in location mean reg model
        # max_iter <- sample(c(2, 10, 100), 1)
        # rel_tol <- runif(1, 1e-6, 1e-5)
        # checking G matrix from cpp and R
        max_iter <- 100
        rel_tol <- 1e-5
        G.cpp <- mr.evalG(y, X, beta0)
        G.R <- mr.evalG_R(y, X, beta0)
        expect_equal(G.R, G.cpp)
        omegahat.cpp <- omega.hat(G = G.cpp, max_iter = max_iter, rel_tol = rel_tol, verbose = FALSE)
        omegahat.R <- omega.hat_R(G = G.R, max_iter = max_iter, rel_tol = rel_tol, verbose = FALSE)
        expect_equal(omegahat.R, omegahat.cpp)
        
        # TODO: may should choose a random percentage of censored observations
        # c <- abs(rnorm(n, mean=y, sd=1)) # independent random censoring values
        # deltas <- y <= c # delta == 1 if observed 0 if censored
        # y[!deltas] <- c[!deltas] # observed lifetime after censoring
        # 1-sum(deltas)/n
        # debug:
        epsilons <- c(y - X %*% beta0)
        deltas <- rep(1,n)
        # censind <- sample(n,2)
        censind <- order(epsilons)[c(4,7)]
        deltas[censind] <- 0
        y[censind] <- rnorm(1,mean=y[censind]) # makes only 1 censored 
        G.cpp <- mr.evalG(y, X, beta0)
        omegahat.cpp <- omega.hat(G.cpp, deltas, epsilons, max_iter = max_iter, rel_tol = rel_tol, verbose = FALSE)
        omegahat.R <- omega.hat_R(G.R, deltas, epsilons, max_iter = max_iter, rel_tol = rel_tol, verbose = FALSE)
        })
        omegahat.cpp
        omegahat.R
        omegahat.cpp - omegahat.R
        # TODO: check the correctness of the R version by the old 5 variables way (?)
        # and then check the c++ version against it 
    }
})
