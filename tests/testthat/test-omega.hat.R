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
        n <- sample(10:20,1)
        p <- sample(1:(n-2), 1)
        X <- replicate(p, rnorm(n))
        beta0 <- rnorm(p)
        y <- c(X %*% beta0) + rnorm(n) # with N(0,1) error term
        lambda0 <- rnorm(p) # m == p in location mean reg model
        max_iter <- sample(c(2, 10, 100), 1)
        rel_tol <- runif(1, 1e-6, 1e-5)
        # checking G matrix from cpp and R
        G.cpp <- mr.evalG(y, X, beta0)
        G.R <- mr.evalG_R(y, X, beta0)
        expect_equal(G.R, G.cpp)
        omegahat.cpp <- omega.hat(G.cpp, max_iter = max_iter, rel_tol = rel_tol, verbose = FALSE)
        omegahat.ROut <- omega.hat_R(G.R, max_iter = max_iter, rel_tol = rel_tol, verbose = FALSE)
        omegahat.R <- omegahat.ROut$omegas
        omegahat.R
        expect_equal(omegahat.R, omegahat.cpp)
        
        # TODO: may should choose a random percentage of censored observations
        c <- abs(rnorm(n, mean=y, sd=1)) # independent random censoring values
        deltas <- y <= c # delta == 1 if observed 0 if censored
        sum(deltas)
        y[!deltas] <- c[!deltas] # observed lifetime after censoring
        epsilons <- c(y - X %*% beta0)
        G.cpp <- mr.evalG(y, X, beta0)
        omegahat.cpp <- omega.hat(G.cpp, epsilons, deltas, max_iter = max_iter, rel_tol = rel_tol, verbose = FALSE)
    }
})
