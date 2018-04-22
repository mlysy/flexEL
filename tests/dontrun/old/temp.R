library(bayesEL)
library(testthat)
source("../testthat/el-utils.R")

n <- 10
mu0 <- 1
# X <- matrix(rep(1,n), n, 1) # each row of X is one observation
X <- matrix(rnorm(n), n, 1)
eps <- rnorm(n) # N(0,1) error term
y <- c(X * mu0) + eps
betaInit <- c(lm(y ~ 1)$coefficients) # TODO: does not work ???
G1 <- mr.evalG(y,X,betaInit)
G2 <- mr.evalG(y,X,round(betaInit,6))
# save(G1, file = "G1.RData")
lambdaNR_R(G1,verbose = T)
lambdaNR_R(G2,verbose = T)
rep(1,nrow(G1)) %*% G1
rep(1,nrow(G2)) %*% G2
t(G1)
t(G2)
expect_equal(G1,G2)
