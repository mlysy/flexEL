# ---- testing lambdaNR ----
library(bayesEL) # always load the package (with library)
library(optimCheck)
# source("el-utils.R")
source("~/bayesEL/tests/testthat/el-utils.R")
# source("~/bayesEL/tests/testthat/mle-check.R")

# library(testthat) # not loaded automatically
context("lambdaNR")

# # G function for location mean regression (old)
# Gmean <- function(y, X, beta) {
#     z <- c(y - X %*% beta)
#     cbind(z, z^2 - 1)
# }

ntest <- 50

# Non-censored case: 
# checking R and C++ implementations are equal, and optimality of the solution 
# TODO: run the for loop ok, but not test_that 
test_that("lambda.R == lambda.cpp", {
    for(ii in 1:ntest) {
        # Location model + mean regression 
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
        
        # Check optimality by optimCheck if converged 
        if (nrout$convergence) {
            # TODO: it seems to be very unstable with p == 1, use optim_proj much better 
            if (p == 1) {
                # xfit <- optim(par = lambda.cpp, fn = Qfun,
                #               # lower = abs(lambda.cpp)*(-1.5), upper = abs(lambda.cpp)*1.5, # with "Brent"
                #               method = "BFGS")
                # ocheck <- optim_refit(xsol = lambda.cpp, Qfun, xopt = xfit$par)
                ocheck <- optim_proj(xsol = lambda.cpp, fun = Qfun, plot = FALSE)
            }
            else {
                ocheck <- optim_refit(xsol = lambda.cpp, Qfun)
            }
            if (max.xdiff(ocheck) > 0.01) print(orefit)
            expect_lt(max.xdiff(ocheck),0.01)
        }
    }
})

# Censored case: 
# TODO: with censoring sometimes only one converges :(
test_that("lambdaC.R == lambdaC.cpp", {
    for(ii in 1:ntest) {
        n <- sample(10:20,1)
        p <- sample(1:(n-2), 1)
        G <- matrix(rnorm(n*p),n,p) # randomly generated G
        max_iter <- sample(c(2, 10, 100), 1)
        rel_tol <- runif(1, 1e-6, 1e-5)
        weights <- abs(rnorm(n))
        weights <- weights / sum(weights) * n # sum(weights) == n
        lambda.cpp <- lambdaNR(G = G, weights,
                               max_iter = max_iter, rel_tol = rel_tol, verbose = FALSE)
        # lambda.cpp
        # ocheck.cpp <- optim_refit(lambda.cpp, QfunCens)
        # expect_lt(max.xdiff(ocheck.cpp),0.01)
        
        nrout <- lambdaNRC_R(G = G, weights,
                             max_iter = max_iter, rel_tol = rel_tol, verbose = FALSE)
        lambda.R <- nrout$lambda
        # lambda.R 
        # ocheck.R <- optim_proj(lambda.R, QfunCens)
        # expect_lt(max.xdiff(ocheck.R),0.01)
        expect_equal(lambda.cpp, lambda.R)
        expect_lt(max.xdiff(ocheck.cpp), 0.01)
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
