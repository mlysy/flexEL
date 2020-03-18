# test evalWeights is working properly
# library(bayesEL) # always load the package (with library)
source("el-utils.R")
source("el-rfuns.R")
source("el-model.R")

# library(testthat) # not loaded automatically
context("evalWeights")

ntest <- 50
test_that("weights.R == weights.cpp", {
  for(ii in 1:ntest) {
    nObs <- sample(10:50,1)
    nEqs <- sample(2:5,1)
    numcens <- sample(1:floor(nObs/2),1)
    censinds <- sample(1:nObs,numcens)
    deltas <- rep(1,nObs)
    deltas[censinds] <- 0
    omegas <- abs(rnorm(nObs))
    omegas <- omegas / sum(omegas) # prob vector
    epsilons <- rnorm(nObs)
    support <- FALSE
    weights.cpp <- .EvalWeights(deltas, omegas, epsilons, support)
    expect_equal(sum(weights.cpp),nObs)
    weights.R <- evalWeights_R(deltas, omegas, epsilons)
    expect_equal(sum(weights.R),nObs)
    expect_equal(weights.cpp, weights.R)
  }
})

ntest <- 50
test_that("weights.R == weights.cpp (with support correction)", {
  for(ii in 1:ntest) {
    nObs <- sample(10:50,1)
    nEqs <- sample(2:5,1)
    numcens <- sample(1:floor(nObs/2),1)
    censinds <- sample(1:nObs,numcens)
    deltas <- rep(1,nObs)
    deltas[censinds] <- 0
    omegas <- abs(rnorm(nObs+1))
    omegas <- omegas / sum(omegas) # prob vector
    epsilons <- rnorm(nObs)
    support <- TRUE
    weights.cpp <- .EvalWeights(deltas, omegas, epsilons, support)
    expect_equal(sum(weights.cpp),nObs+1)
    weights.R <- evalWeights_R(c(deltas,0), omegas, c(epsilons,-Inf))
    expect_equal(sum(weights.R),nObs+1)
    expect_equal(weights.cpp, weights.R)
  }
})

