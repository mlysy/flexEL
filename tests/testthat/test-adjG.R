source("../testthat/el-utils.R")

# library(testthat) # not loaded automatically
context("adjG")

ntest <- 50

test_that("adjG.R == adjR.cpp", {
  for(ii in 1:ntest) {
    n <- sample(1:100,1)
    p <- sample(1:100,1)
    max_iter <- 500
    rel_tol <- 1e-7
    G <- matrix(rnorm(n*p), n, p)
    a <- abs(rnorm(1))
    aG.R <- adjG_R(G,a)
    aG.cpp <- adjG(G,a)
    expect_equal(aG.R, aG.cpp)
  }
})
