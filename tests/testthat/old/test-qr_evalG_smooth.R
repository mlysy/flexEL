#---- testing qr_evalG (location model) ----

source("el_rfuns.R")
source("reg_models.R")

# library(testthat) # not loaded automatically
context("qr_evalG_smooth")

ntest <- 5

test_that("qr_evalG_R == qr_evalG_cpp (one quantile level)", {
  for(ii in 1:ntest) {
    n <- sample(10:20,1)
    p <- sample(1:(n-2), 1)
    X <- replicate(p, rnorm(n))
    X[,1] <- rep(1,n)
    beta <- rnorm(p)
    alpha <- runif(1)
    y <- c(X %*% beta) + rnorm(n) # with N(0,1) error term
    sp <- sample(1:100, 1)
    # checking G matrix from cpp and R
    G_cpp <- flexEL::qr_evalG(y, X, alpha, beta, s = sp)
    G_R <- qr_evalG_smooth_R(y, X, alpha, beta, s = sp)
    expect_equal(G_cpp,G_R)
  }
})

reorg_jacobian <- function(jac, n, p) {
  ll <- split(split(jac, row(jac)), rep(seq(from = 1, to = n), p))
  lm <- lapply(ll, function(l) {
    matrix(unlist(l), nrow = p)
  })
  do.call(rbind, lm)
}

test_that("smooth qr dGdt matches numDeriv's jacobian (one quantile level)", {
  for(ii in 1:ntest) {
    n <- sample(10:20,1)
    p <- sample(1:(n-2), 1)
    X <- replicate(p, rnorm(n))
    X[,1] <- rep(1,n)
    beta <- rnorm(p)
    alpha <- runif(1)
    y <- c(X %*% beta) + rnorm(n) # with N(0,1) error term
    sp <- sample(1:100, 1)
    dGdt_cpp <- flexEL:::QuantReg_dGdt_smooth(y, t(X), c(1,alpha), matrix(beta, ncol = 1), sp)
    qr_evalG_fix <- function(bb) {
      flexEL::qr_evalG(y, X, alpha, bb, sp)
    }
    dGdt_R <- numDeriv::jacobian(qr_evalG_fix, beta, method.args = list(eps = 1e-5))
    expect_equal(dGdt_cpp, reorg_jacobian(dGdt_R, n, p))
  }
})
