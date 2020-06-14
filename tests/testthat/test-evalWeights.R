# ---- test evalWeights is working properly ----

source("el_rfuns.R")

# library(testthat) # not loaded automatically
context("evalWeights")

ntest <- 50
test_that("weights_R == weights_cpp", {
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
    weights_cpp <- flexEL:::.EvalWeights(deltas, omegas, epsilons, support)
    expect_equal(sum(weights_cpp),nObs)
    weights_R <- evalWeights_R(deltas, omegas, epsilons)
    expect_equal(sum(weights_R),nObs)
    expect_equal(weights_cpp, weights_R)
  }
})

ntest <- 50
test_that("weights_R == weights_cpp (with support correction)", {
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
    weights_cpp <- flexEL:::.EvalWeights(deltas, omegas, epsilons, support)
    expect_equal(sum(weights_cpp),nObs+1)
    weights_R <- evalWeights_R(c(deltas,0), omegas, c(epsilons,-Inf))
    expect_equal(sum(weights_R),nObs+1)
    expect_equal(weights_cpp, weights_R)
  }
})

