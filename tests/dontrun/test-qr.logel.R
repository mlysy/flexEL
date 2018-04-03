#---- testing qr.logel  ----
library(bayesEL) # always load the package (with library)
# source("el-utils.R")
source("~/bayesEL/tests/testthat/el-utils.R")

# library(testthat) # not loaded automatically
context("qr.logel")

ntest <- 50

test_that("qr.logel.R == qr.logel.cpp", {
    for(ii in 1:ntest) {
        n <- sample(10:20,1)
        p <- sample(1:(n-2), 1)
        X <- replicate(p, rnorm(n))
        X[1,] <- rep(1,p)
        beta <- rnorm(p)
        alpha <- runif(1)
        y <- c(X %*% beta) + rnorm(n) # with N(0,1) error term
        max_iter <- sample(c(2, 10, 100), 1)
        rel_tol <- runif(1, 1e-6, 1e-5)
        # checking logel matrix from cpp and R
        logel.cpp <- qr.logel(y, X, alpha, beta, max_iter = max_iter, rel_tol = rel_tol)
        logel.R <- qr.logel_R(y, X, alpha, beta, max_iter = max_iter, rel_tol = rel_tol)
        expect_equal(logel.cpp,logel.R)
    }
})
