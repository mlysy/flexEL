# ---- testing R and C++ implementations of evalWeights_smooth are equal ----

source("el_rfuns.R")

# library(testthat) # not loaded automatically
context("evalWeights_smooth")

ntest <- 50

test_that("evalPsos_R == evalPsos_cpp", {
  for(ii in 1:ntest) {
    n <- sample(1:100,1)
    s <- sample(1:100,1)
    omegas <- abs(rnorm(n))
    omegas <- omegas/sum(omegas)
    epsilons <- rnorm(n)
    deltas <- rep(1,n)
    numcens <- sample(1:round(n/2),1)
    deltas[sample(1:n,numcens)] <- 0
    support <- FALSE
    val_cpp <- flexEL:::.EvalWeightsSmooth(deltas,omegas,epsilons,s,support)
    val_R <- evalWeights_smooth_R(deltas,omegas,epsilons,s)
    expect_equal(val_R, val_cpp)
  }
})

test_that("evalPsos_R == evalPsos_cpp (with support correction)", {
  for(ii in 1:ntest) {
    n <- sample(1:100,1)
    s <- sample(1:100,1)
    omegas <- abs(rnorm(n+1))
    omegas <- omegas/sum(omegas)
    epsilons <- rnorm(n)
    deltas <- rep(1,n)
    numcens <- sample(1:round(n/2),1)
    deltas[sample(1:n,numcens)] <- 0
    support <- TRUE
    val_cpp <- flexEL:::.EvalWeightsSmooth(deltas,omegas,epsilons,s,support)
    val_R <- evalWeights_smooth_R(c(deltas,0), omegas, c(epsilons,-Inf),s,support)
    expect_equal(val_R, val_cpp)
  }
})
