# ---- testing qrls.evalG implementations in R and C++ are equal ----

source("el_funs.R")
source("reg_models.R")

context("qrls.evalG")

ntest <- 50

test_that("qrls.evalG.R == qrls.evalG.cpp", {
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
    G.cpp <- qrls_evalG(y,X,Z,alpha,beta,gamma,sig2,nu)
    G.R <- qrls_evalG_R(y,X,Z,alpha,beta,gamma,sig2,nu)
    expect_equal(G.cpp, G.R)
  }
})

