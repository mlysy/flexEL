#---- testing qr.eval.G  ----
library(bayesEL) # always load the package (with library)
# source("el-utils.R")
source("~/bayesEL/tests/testthat/el-utils.R")

# library(testthat) # not loaded automatically
context("qr.evalG")

ntest <- 50

test_that("qr.evalG.R == qr.evalG.cpp", {
    for(ii in 1:ntest) {
        n <- sample(10:20,1)
        p <- sample(1:(n-2), 1)
        X <- replicate(p, rnorm(n))
        X[1,] <- rep(1,p)
        beta <- rnorm(p)
        alpha <- runif(1)
        y <- c(X %*% beta) + rnorm(n) # with N(0,1) error term
        # checking G matrix from cpp and R
        G.cpp <- qr.evalG(y, X, alpha, beta)
        G.R <- qr.evalG_R(y, X, alpha, beta)
        expect_equal(G.cpp,G.R)
    }
})
