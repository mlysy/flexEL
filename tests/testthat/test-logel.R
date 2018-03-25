#---- testing logel  ----
library(bayesEL) # always load the package (with library)
# library(optimCheck)
# source("el-utils.R")
source("~/bayesEL/tests/testthat/el-utils.R")

# library(testthat) # not loaded automatically
context("logel")

ntest <- 50

# Non-censored case: 
# checking R and C++ implementations are equal
test_that("logel.R == logel.cpp", {
    for(ii in 1:ntest) {
        n <- sample(10:20,1)
        p <- sample(1:(n-2), 1)
        max_iter <- 100
        # max_iter <- sample(c(2, 10, 100), 1)
        rel_tol <- runif(1, 1e-6, 1e-5)
        G <- matrix(rnorm(n*p),n,p) 
        logEL(G, max_iter = max_iter, rel_tol = rel_tol, verbose = TRUE)
        
    }
})