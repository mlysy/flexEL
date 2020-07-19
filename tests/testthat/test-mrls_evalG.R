# ---- testing qrls_evalG implementations in R and C++ are equal ----

source("el_rfuns.R")
source("reg_models.R")

context("mrls_evalG")

ntest <- 50

test_that("mrls_evalG_R == mrls_evalG_cpp", {
  for(ii in 1:ntest) {
    # Location-scale model + mean regression
    n <- sample(10:20,1)
    p <- sample(1:(n-2), 1)
    q <- sample(1:(n-2), 1)
    X <- replicate(p, rnorm(n))
    Z <- replicate(q, rnorm(n))
    beta <- rnorm(p)
    gamma <- rnorm(q)
    sig2 <- abs(rnorm(1)) # scale parameter
    y <- c(X %*% beta + exp(Z %*% gamma)) + rnorm(n) # with N(0,1) error term
    max_iter <- sample(c(2, 10, 100), 1)
    rel_tol <- runif(1, 1e-6, 1e-5)
    # checking G matrix from cpp and R
    G_cpp <- flexEL::mrls_evalG(y,X,Z,beta,gamma,sig2)
    G_R <- mrls_evalG_R(y,X,Z,beta,gamma,sig2)
    expect_equal(G_R, G_cpp)
  }
})

