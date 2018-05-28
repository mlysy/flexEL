# test that hlm is working properly
source("../../R/hlm.R")
source("el-utils.R")
library(optimCheck)
library(numDeriv)

# library(testthat) # not loaded automatically
context("hlm")

ntest <- 50

test_that("hlm optimality", {
  for (ii in 1:ntest) {
   n <- sample(10:50,1)
   p <- sample(1:round(n/2), 1)
   # p <- 3
   q <- sample(1:round(n/2), 1)
   # q <- 2
   X <- matrix(rnorm(n*p),n,p)
   W <- matrix(rnorm(n*q),n,q)
   beta0 <- rnorm(p)
   gamma0 <- rnorm(q)/2
   eps <- rnorm(n)
   yy <- X %*% beta0 + exp(0.5*W %*% gamma0)*eps
   cc <- rnorm(n, mean=2)
   delta <- yy <= cc
   sum(1-delta)/length(delta)
   y <- yy
   y[!delta] <- cc[!delta]
   hlmout <- hlm(y,delta,X,W,rel_tol = 1e-5)
   if (hlmout$conv) {
     ocheck_beta <- optim_proj(xsol = hlmout$coef$beta,
                               fun = function(beta) hlm_loglik(beta,
                                                               hlmout$coef$gamma,
                                                               y,delta,X,W),
                               plot = FALSE)
     expect_lt(max.xdiff(ocheck_beta),0.01)
     ocheck_gamma <- optim_proj(xsol = hlmout$coef$gamma,
                                fun = function(gamma) hlm_loglik(hlmout$coef$beta,
                                                                 gamma,
                                                                 y,delta,X,W),
                                plot = FALSE)
     expect_lt(max.xdiff(ocheck_gamma),0.01)
   }
  }
})
