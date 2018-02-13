library(bayesEL) # always load the package (with library)
source("el-utils.R")
## source("../bayesEL/tests/testthat/el-utils.R")
## source("../bayesEL/tests/donotrun/mle-check.R")

# library(testthat) # not loaded automatically
context("lambdaNR")

# G function for mean regression
Gmean <- function(y, X, beta) {
    z <- c(y - X %*% beta)
    cbind(z, z^2 - 1)
}

ntest <- 50

test_that("lambda.R == lambda.cpp", {
    for(ii in 1:ntest) {
        #  simulate data
        m <- 2 # for this formulation always have m = 2
        n <- sample(10:20,1)
        p <- sample(1:(n-2), 1)
        X <- replicate(p, rnorm(n))
        beta0 <- rnorm(p)
        y <- c(X %*% beta0) + rnorm(n) # with N(0,1) error term
        lambda0 <- rnorm(m) #
        max_iter <- sample(c(2, 10, 100), 1)
        eps <- runif(1, 1e-6, 1e-5)
        G <- Gmean(y, X, beta0)
        lambda.cpp <- lambdaNR(G = G,
                               max_iter = max_iter, eps = eps, verbose = FALSE)
        # lambda.cpp
        # now in R
        nrout <- lambdaNR_R(G = G,
                            max_iter = max_iter, eps = eps)
        lambda.R <- nrout$lambda
        expect_equal(lambda.R, lambda.cpp)
    }
})

#-----------------------------------------------------------------------------------------------------------------------------------------------------------

## # some of these inputs can be processed in R, not C++
## # usually each row of X is an obs, so have R code transpose it for C++
## ll.cpp <- logELmean(nObs = n, nEqs = 2,
##                     y = y, X = t(X),
##                     beta = beta0, lambda0 = beta0)
