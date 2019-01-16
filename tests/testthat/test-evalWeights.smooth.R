# ---- testing R and C++ implementations of evalWeights.smooth are equal ----

source("../dontrun/smoothEL.R")
source("../dontrun/mode-functions.R")

# library(testthat) # not loaded automatically
context("evalWeights.smooth")

ntest <- 50

test_that("evalPsos.R == evalPsos.cpp", {
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
    val.cpp <- .evalWeights.smooth(deltas,omegas,epsilons,s,support)
    val.R <- evalWeights.smooth_R(deltas,omegas,epsilons,s)
    expect_equal(val.R, val.cpp)
  }
})
