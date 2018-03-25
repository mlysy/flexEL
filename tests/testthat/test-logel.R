#---- testing logel  ----
library(bayesEL) # always load the package (with library)
# library(optimCheck)
# source("el-utils.R")
source("~/bayesEL/tests/testthat/el-utils.R")

# library(testthat) # not loaded automatically
context("logEL")

ntest <- 50

# Non-censored case: 
# checking R and C++ implementations are equal 
test_that("logel.R == logel.cpp", {
    for(ii in 1:ntest) {
        n <- sample(10:20,1)
        p <- sample(1:(n-5), 1)
        max_iter <- 100
        # max_iter <- sample(c(2, 10, 100), 1)
        # rel_tol <- runif(1, 1e-6, 1e-5)
        rel_tol <- 1e-6
        G <- matrix(rnorm(n*p),n,p)
        omegas <- omega.hat(G) # optimal omegas 
        if (sum(omegas) != 0) {
            logopt.cpp <- logEL(omegas, G)
            logopt.R <- logEL_R(omegas, G)
            expect_equal(logopt.cpp,logopt.R)
        }
        omegas <- omegas + rnorm(n,sd=0.01)
        omegas <- abs(omegas) / sum(abs(omegas))
        if (sum(omegas) != 0) {
            logtwk.cpp <- logEL(omegas, G)
            logtwk.R <- logEL_R(omegas, G)
            expect_equal(logtwk.cpp,logtwk.R)
        }
        expect_gt(logopt - logtwk, -0.01) # TODO: ??
        # X <- matrix(rnorm(n*p),n,p)
        # beta <- 2+ rnorm(p)
        # y <- c(X %*% beta + rnorm(n))
        # G <- mr.evalG(y,X,beta)
        # Old code: find the maximized loglikelihood
        # logEL(G, max_iter = max_iter, rel_tol = rel_tol, verbose = TRUE)
    }
})
