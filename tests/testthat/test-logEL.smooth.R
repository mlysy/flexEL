# ---- testing R and C++ implementations of logEL.smooth are equal ----

source("../dontrun/smoothEL.R")

# library(testthat) # not loaded automatically
context("logEL.smooth")

ntest <- 50

test_that("logEL.smooth.R == logEL.smooth.cpp", {
  for(ii in 1:ntest) {
    n <- sample(1:100,1)
    s <- sample(1:100,1)
    omegas <- abs(rnorm(n))
    omegas <- omegas/sum(omegas)
    epsilons <- rnorm(n)
    deltas <- rep(1,n)
    numcens <- sample(1:round(n/2),1)
    deltas[sample(1:n,numcens)] <- 0
    val.cpp <- logEL.smooth(omegas,epsilons,deltas,s)
    val.R <- logEL.smooth_R(omegas,epsilons,deltas,s)
    expect_equal(val.R, val.cpp)
  }
})
