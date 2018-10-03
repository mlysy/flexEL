# ---- testing qrls.evalG implementations in R and C++ are equal ----

# library(bayesEL)
source("el-utils.R")
source("el-rfuns.R")
source("el-model.R")
source("../dontrun/mode-functions.R")
# library(testthat)

context("mrls.evalG")

ntest <- 50

test_that("qrls.evalG.smooth.R == qrls.evalG.smooth.cpp", {
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
    nu <- rnorm(1)
    sig2 <- abs(rnorm(1))
    # sig2 <- 1
    y <- c(X %*% beta + sqrt(sig2)*exp(Z %*% gamma)*rnorm(n)) # with multiplicative N(0,1) error
    # checking G matrix from cpp and R
    sp <- sample(100,1)
    G.cpp <- qrls.evalG.smooth(y,X,Z,alpha,beta,gamma,sig2,nu,sp)
    G.R <- qrls.evalG.smooth_R(y,X,Z,alpha,beta,gamma,sig2,nu,sp)
    expect_equal(G.cpp, G.R)
  }
})

