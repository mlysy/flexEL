#---- testing logel  ----
library(bayesEL) # always load the package (with library)
# library(optimCheck)
# source("el-utils.R")
source("~/bayesEL/tests/testthat/el-utils.R")

# library(testthat) # not loaded automatically
context("logEL")

ntest <- 50

# Non-censored case: 
# checking R and C++ implementations are equal
# check the optimal omegas obtained is better than a random feasible omegas_twk 
#   which is chosen from the null space of G.
test_that("logel.R == logel.cpp", {
    for(ii in 1:ntest) {
        n <- sample(10:20,1)
        p <- sample(1:(n-5), 1)
        max_iter <- 100
        # max_iter <- sample(c(2, 10, 100), 1)
        # rel_tol <- runif(1, 1e-6, 1e-5)
        rel_tol <- 1e-6
        G <- matrix(rnorm(n*p),n,p)
        omegas <- omega.hat(G) # optimal omegas (by C++ here)
        
        if (sum(omegas) != 0) {
            # omegas should satisfy: omegas %*% G = 0
            expect_equal(c(omegas %*% G),rep(0,p)) 
            
            # check C++ and R are the same 
            logopt.cpp <- logEL(omegas, G)
            logopt.R <- logEL_R(omegas, G)
            expect_equal(logopt.cpp,logopt.R)
            
            # use Null space of G for checking the optimality of omegas
            NG <- Null(G) # any col c of NG satisfies c %*% G = 0
            # we know that omegas > 0, just take any c of NG, add omegas until > 0
            # then normalize it
            ind <- sample(n-p,1) # dim of NG is n x (n-p)
            omegas_twk <- NG[,ind] + 50*omegas
            # make sure omegas_twk is nonnegative: 
            while (sum(omegas_twk < 0) > 0) omegas_twk <- omegas_twk + 50*omegas
            omegas_twk <- omegas_twk / sum(omegas_twk) # normalize it
            # this should be 0s still 
            expect_equal(c(omegas_twk %*% G), rep(0,p))
            logtwk.cpp <- logEL(omegas_twk, G)
            # logopt should be greater than logtwk 
            expect_gt(logopt - logtwk, 0)
        }
 
        # randomly tweak omegas
        # omegas <- omegas + rnorm(n,sd=0.01)
        # omegas <- abs(omegas) / sum(abs(omegas))
        # if (sum(omegas) != 0) {
        #     logtwk.cpp <- logEL(omegas, G)
        #     logtwk.R <- logEL_R(omegas, G)
        #     expect_equal(logtwk.cpp,logtwk.R)
        # }
        # expect_gt(logopt - logtwk, -0.01) # TODO: ??
        
        # location model data generation 
        # X <- matrix(rnorm(n*p),n,p)
        # beta <- 2+ rnorm(p)
        # y <- c(X %*% beta + rnorm(n))
        # G <- mr.evalG(y,X,beta)
        # Old code: find the maximized loglikelihood
        # logEL(G, max_iter = max_iter, rel_tol = rel_tol, verbose = TRUE)
    }
})
