# ---- testing qrls.evalG implementations in R and C++ are equal ----

source("el-utils.R")
# library(testthat)

context("mrls.evalG")

ntest <- 50

test_that("mr.evalG.R == mr.evalG.cpp", {
  for(ii in 1:ntest) {
    # Location-scale model + quantile regression
    n <- sample(10:20,1)
    p <- sample(1:(n-2), 1)
    q <- sample(1:(n-2), 1)
    X <- replicate(p, rnorm(n))
    Z <- replicate(q, rnorm(n))
    alpha <- runif(1)
    beta <- rnorm(p)
    gamma <- rnorm(q)
    y <- c(X %*% beta + exp(Z %*% gamma)) + rnorm(n) # with N(0,1) error term
    max_iter <- sample(c(2, 10, 100), 1)
    rel_tol <- runif(1, 1e-6, 1e-5)
    # checking G matrix from cpp and R
    G.cpp <- qrls.evalG(y,X,Z,alpha,beta,gamma)
    G.R <- qrls.evalG_R(y,X,Z,alpha,beta,gamma)
    expect_equal(G.cpp, G.R)
  }
})


