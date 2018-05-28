# test evalEpsilonsLS is working properly
# library(bayesEL) # always load the package (with library)
source("el-utils.R")
source("el-rfuns.R")
source("el-model.R")

# library(testthat) # not loaded automatically
context("evalEpsilons")

ntest <- 50
test_that("epsilons.R == epsilons.cpp", {
  for(ii in 1:ntest) {
    nObs <- sample(10:50,1)
    nBet <- sample(2:5,1)
    nGam <- sample(2:5,1)
    beta <- rnorm(nBet)
    gamma <- rnorm(nGam)
    sig2 <- abs(rnorm(1,mean=1))
    X <- matrix(rnorm(nObs*nBet),nObs,nBet)
    Z <- matrix(rnorm(nObs*nGam),nObs,nGam)
    y <- rnorm(nObs)
    epsilons.cpp <- evalEpsilonsLS(y,X,Z,beta,gamma,sig2)
    epsilons.R <- c(evalEpsilonsLS_R(y,X,Z,beta,gamma,sig2))
  }
})
