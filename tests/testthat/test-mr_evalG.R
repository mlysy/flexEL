#---- testing mr_evalG (location model) ----

source("reg_models.R")

# library(testthat) # not loaded automatically
context("mr_evalG")

ntest <- 5

test_that("mr_evalG_R == mr_evalG_cpp", {
  for(ii in 1:ntest) {
    # Location model + mean regression
    n <- sample(10:20,1)
    p <- sample(1:(n-2), 1)
    X <- replicate(p, rnorm(n))
    beta0 <- rnorm(p)
    y <- c(X %*% beta0) + rnorm(n) # with N(0,1) error term
    # checking G matrix from cpp and R
    G_cpp <- flexEL::mr_evalG(y, X, beta0)
    G_R <- mr_evalG_R(y, X, beta0)
    expect_equal(G_R, G_cpp)
  }
})

reorg_jacobian <- function(jac, n, p) {
  ll <- split(split(jac, row(jac)), rep(seq(from = 1, to = n), p))
  lm <- lapply(ll, function(l) {
    matrix(unlist(l), nrow = p)
  })
  do.call(rbind, lm)
}

test_that("mr dGdt matches numDeriv's jacobian", {
  for(ii in 1:ntest) {
    # Location model + mean regression
    n <- sample(10:20,1)
    p <- sample(1:(n-2), 1)
    X <- replicate(p, rnorm(n))
    beta0 <- rnorm(p)
    y <- c(X %*% beta0) + rnorm(n) # with N(0,1) error term
    dGdt_cpp <- flexEL:::.MeanRegEvaldGdt(y, t(X), beta0)
    mr_evalG_fix <- function(bb) {
      flexEL::mr_evalG(y, X, bb)
    }
    dGdt_R <- numDeriv::jacobian(mr_evalG_fix, beta0, method.args = list(eps = 1e-5))
    expect_equal(dGdt_cpp, reorg_jacobian(dGdt_R, n, p))
  }
})
