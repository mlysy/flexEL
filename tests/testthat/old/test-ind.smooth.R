# ---- testing R and C++ implementations of omega.hat are equal and optimality ----
# library(bayesEL) # always load the package (with library)
# source("el-utils.R")
# source("el-rfuns.R")
# source("el-model.R")
source("../dontrun/smoothEL.R")

# library(testthat) # not loaded automatically
context("ind.smooth")

ntest <- 50

test_that("ind.smooth.R == ind.smooth.cpp", {
  for(ii in 1:ntest) {
    x <- rnorm(1)
    s <- sample(1:100,1)
    val.R <- ind.smooth_R(x,s)
    val.cpp <- .ind.smooth(x,s)
    expect_equal(val.R, val.cpp)
  }
})

test_that("ind1.smooth.R == ind1.smooth.cpp", {
  for(ii in 1:ntest) {
    x <- rnorm(1)
    s <- sample(1:100,1)
    val.R <- ind1.smooth_R(x,s)
    val.cpp <- .ind1.smooth(x,s)
    expect_equal(val.R, val.cpp)
  }
})
