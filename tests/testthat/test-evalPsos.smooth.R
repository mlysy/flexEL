# ---- testing R and C++ implementations of evalPsos are equal ----

source("../dontrun/smoothEL.R")

# library(testthat) # not loaded automatically
context("evalPsos")

ntest <- 50

test_that("evalPsos.R == evalPsos.cpp", {
  for(ii in 1:ntest) {
    n <- sample(1:20,1)
    s <- sample(1:100,1)
    omegas <- abs(rnorm(n))
    omegas <- omegas/sum(omegas)
    epsilons <- rnorm(n)
    ii <- sample(1:n,1)
    val.cpp <- evalPsos.smooth(ii,omegas,epsilons,s)
    val.R <- evalPsos.smooth_R(ii,omegas,epsilons,s)
    expect_equal(val.R, val.cpp)
  }
})
