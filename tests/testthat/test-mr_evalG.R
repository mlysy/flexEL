#---- testing mr.evalG (location model) ----
## library(bayesEL) # always load the package (with library)
source("el-utils.R")
source("el-rfuns.R")
source("el-model.R")
## source("~/bayesEL/tests/testthat/el-utils.R")

# library(testthat) # not loaded automatically
context("mr_evalG")

ntest <- 50

test_that("mr.evalG.R == mr.evalG.cpp", {
  for(ii in 1:ntest) {
    # Location model + mean regression
    n <- sample(10:20,1)
    p <- sample(1:(n-2), 1)
    X <- replicate(p, rnorm(n))
    beta0 <- rnorm(p)
    y <- c(X %*% beta0) + rnorm(n) # with N(0,1) error term
    # checking G matrix from cpp and R
    G_cpp <- mr_evalG(y, X, beta0)
    G_R <- mr.evalG_R(y, X, beta0)
    expect_equal(G_R, G_cpp)
  }
})

