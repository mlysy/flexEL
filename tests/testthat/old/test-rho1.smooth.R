# ---- testing R and C++ implementations of rho1.smooth are equal ----

source("../dontrun/smoothEL.R")

# library(testthat) # not loaded automatically
context("evalWeights.smooth")

ntest <- 50

test_that("evalPsos.R == evalPsos.cpp", {
  for(ii in 1:ntest) {
    s <- sample(1:100,1)
    x <- rnorm(1)
    tau <- runif(1)
    val.cpp <- .rho1.smooth(x,tau,s)
    val.R <- rho1.smooth_R(x,tau,s)
    expect_equal(val.R, val.cpp)
  }
})
