#---- testing qr_evalG (location model) ----

source("el_rfuns.R")
source("reg_models.R")

# library(testthat) # not loaded automatically
context("qr_evalG")

ntest <- 50

test_that("qr_evalG_R == qr_evalG_cpp", {
    for(ii in 1:ntest) {
        n <- sample(10:20,1)
        p <- sample(1:(n-2), 1)
        X <- replicate(p, rnorm(n))
        X[1,] <- rep(1,p)
        beta <- rnorm(p)
        alpha <- runif(1)
        y <- c(X %*% beta) + rnorm(n) # with N(0,1) error term
        # checking G matrix from cpp and R
        G_cpp <- qr_evalG(y, X, alpha, beta)
        G_R <- qr_evalG_R(y, X, alpha, beta)
        expect_equal(G_cpp,G_R)
    }
})
