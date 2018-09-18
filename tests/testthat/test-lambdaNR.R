# ---- testing lambdaNR ----
# library(bayesEL) # always load the package (with library)
library(optimCheck)
source("el-utils.R")
source("el-rfuns.R")
source("el-model.R")

# library(testthat) # not loaded automatically
context("lambdaNR")

ntest <- 50

# Non-censored case:
# checking R and C++ implementations are equal
test_that("no censoring: lambda.R == lambda.cpp", {
  for(ii in 1:ntest) {
    n <- sample(10:20,1)
    p <- sample(1:(n-2), 1)
    G <- matrix(rnorm(n*p),n,p) # randomly generated G
    max_iter <- sample(c(2, 10, 100), 1)
    rel_tol <- runif(1, 1e-6, 1e-5)
    lambda.cpp <- lambdaNR(G = G,
                           max_iter = max_iter, rel_tol = rel_tol, verbose = FALSE)
    nrout <- lambdaNR_R(G = G,
                        max_iter = max_iter, rel_tol = rel_tol)
    lambda.R <- nrout$lambda
    expect_equal(lambda.R, lambda.cpp)
  }
})

# checking optimality of the solution from C++
test_that("no censoring: lambda.cpp is optimal", {
  for(ii in 1:ntest) {
    n <- sample(10:20,1)
    p <- sample(1:(n-2), 1)
    G <- matrix(rnorm(n*p),n,p) # randomly generated G
    max_iter <- sample(c(5, 10, 100), 1)
    rel_tol <- runif(1, 1e-6, 1e-5)
    lambda.cpp <- lambdaNR(G = G,
                           max_iter = max_iter, rel_tol = rel_tol, verbose = FALSE)
    # check optimality by optim_proj if converged
    if (!any(is.na(lambda.cpp))) {
      ocheck <- optim_proj(xsol = lambda.cpp,
                           fun = function(lambda) Qfun(lambda, G),
                           plot = FALSE)
      expect_lt(max.xdiff(ocheck),0.01)
    }
  }
})

# lambdahat <- lambda.R
# G <- G.R
# par(mfrow=c(1,1))
# curve(sapply(x, QfunCens, G = G), from = lambdahat-1, to = lambdahat+1)
# abline(v = lambdahat, col='red')

#-----------------------------------------------------------------------------------------------------------------------------------------------------------

## # some of these inputs can be processed in R, not C++
## # usually each row of X is an obs, so have R code transpose it for C++
## ll.cpp <- logELmean(nObs = n, nEqs = 2,
##                     y = y, X = t(X),
##                     beta = beta0, lambda0 = beta0)
