# ---- testing logel implementation in R and C++ are equal ----

source("el_rfuns.R")

# library(testthat) # not loaded automatically
context("logEL")

ntest <- 50

# Non-censored case:
test_that("logel_R == logel_cpp", {
  for(ii in 1:ntest) {
    n <- sample(10:20,1)
    omegas <- abs(rnorm(n))
    omegas <- omegas/sum(omegas)
    logopt_cpp <- logEL(omegas)
    logopt_R <- logEL_R(omegas)
    expect_equal(logopt_cpp,logopt_R)
  }
})

# Censored case:
test_that("logelC_R == logelC_cpp", {
  for(ii in 1:ntest) {
    n <- sample(10:20,1)
    omegas <- abs(rnorm(n))
    k <- sample(round(n/2),1)
    idx <- sample(n,k)
    omegas[idx] <- omegas[idx]/(10^sample(5:25,k))
    omegas <- omegas/sum(omegas)
    deltas <- rep(1,n)
    deltas[idx] <- 0
    epsilons <- rnorm(n)
    logopt_cpp <- logEL(omegas,epsilons,deltas)
    logopt_R <- logEL_R(omegas,epsilons,deltas)
    expect_equal(logopt_cpp,logopt_R)
  }
})
