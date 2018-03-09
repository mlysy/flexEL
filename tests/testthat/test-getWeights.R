# test evalWeights is working properly
library(bayesEL) # always load the package (with library)
source("~/bayesEL/tests/testthat/el-utils.R")

# library(testthat) # not loaded automatically
context("getWeights")

ntest <- 50
test_that("lambda.R == lambda.cpp", {
    for(ii in 1:ntest) {
        nObs <- sample(10:50,1) 
        nEqs <- sample(2:5,1)
        X <- matrix(rnorm(nObs*nEqs),nObs,nEqs)
        beta <- rnorm(nEqs)
        y <- X %*% beta + rnorm(nObs)
        numcens <- sample(1:floor(nObs/2),1)
        censinds <- sample(1:nObs,numcens)
        deltas <- rep(1,nObs)
        deltas[censinds] <- 0
        omegas <- abs(rnorm(nObs)) 
        omegas <- omegas / sum(omegas) # prob vector
        weights.R <- getWeights_R(y, X, deltas, omegas, beta)
        expect_equal(sum(weights.R),nObs)
        weights.cpp <- getWeights(y, X, deltas, omegas, beta)
        expect_equal(sum(weights.cpp),nObs)
    }
})
