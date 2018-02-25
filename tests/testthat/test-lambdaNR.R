#---- testing evalG and lambdaNR (location + mean regression + censoring) ----
library(bayesEL) # always load the package (with library)
# source("el-utils.R")
source("~/bayesEL/tests/testthat/el-utils.R")
source("~/bayesEL/tests/testthat/mle-check.R")

# library(testthat) # not loaded automatically
context("lambdaNR")

# # G function for mean regression (old)
# Gmean <- function(y, X, beta) {
#     z <- c(y - X %*% beta)
#     cbind(z, z^2 - 1)
# }

ntest <- 50

test_that("lambda.R == lambda.cpp", {
    for(ii in 1:ntest) {
        # Location model + mean regression 
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
        # checking lambda from cpp and R
        lambda.cpp <- lambdaNR(G = G.cpp,
                               max_iter = max_iter, rel_tol = rel_tol, verbose = FALSE)

        nrout <- lambdaNR_R(G = G.R,
                            max_iter = max_iter, rel_tol = rel_tol)
        lambda.R <- nrout$lambda
        expect_equal(lambda.R, lambda.cpp)
        
        # TODO: Location model + mean regression + censoring
    }
})

#-----------------------------------------------------------------------------------------------------------------------------------------------------------

## # some of these inputs can be processed in R, not C++
## # usually each row of X is an obs, so have R code transpose it for C++
## ll.cpp <- logELmean(nObs = n, nEqs = 2,
##                     y = y, X = t(X),
##                     beta = beta0, lambda0 = beta0)
