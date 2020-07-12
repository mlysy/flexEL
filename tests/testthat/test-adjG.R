
source("el_rfuns.R")

# library(testthat) # not loaded automatically
context("adjG")

ntest <- 50

test_that("adjG_R == adjR_cpp", {
  for(ii in 1:ntest) {
    n <- sample(1:100,1)
    p <- sample(1:100,1)
    G <- matrix(rnorm(n*p), n, p)
    a <- abs(rnorm(1))
    aG_cpp <- flexEL:::adjG(G,a)
    aG_R <- adjG_R(G,a)
    expect_equal(aG_R, aG_cpp)
  }
})
