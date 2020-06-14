# ---- testing R and C++ implementations of logEL_smooth are equal ----

source("el_rfuns.R")

# library(testthat) # not loaded automatically
context("logEL_smooth")

ntest <- 50

test_that("logEL_smooth_R == logEL_smooth_cpp", {
  for(ii in 1:ntest) {
    n <- sample(1:100,1)
    s <- sample(1:100,1)
    omegas <- abs(rnorm(n))
    omegas <- omegas/sum(omegas)
    epsilons <- rnorm(n)
    deltas <- rep(1,n)
    numcens <- sample(1:round(n/2),1)
    deltas[sample(1:n,numcens)] <- 0
    val_cpp <- logEL_smooth(omegas,epsilons,deltas,s)
    val_R <- logEL_smooth_R(omegas,epsilons,deltas,s)
    expect_equal(val_R, val_cpp)
  }
})
