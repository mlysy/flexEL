#--- check that lambdaNR.R is working properly ---------------------------------

require(bayesELnew)
source("../bayesELnew/tests/testthat/el-utils.R")
source("../bayesELnew/tests/donotrun/mle-check.R")

# 1-d problem
N <- 10 # number of observations
m <- 1 # number of dimensions
y <- rnorm(N*m)
X <- matrix(rep(1,N))
G <- matrix(y, m, N)

# optimization in C++
lambda0 <- rnorm(m)
lambdahat <- lambdaNR(y = y, X = X, G = t(G), nObs = N, nEqs = m, lambda0 = lambda0)

# by implementation in R
lambdahat_Rout <- lambdaNR_R(G = t(G), lambda0 = lambda0)
lambdahat_R <- lambdahat_Rout$lambda

# difference of the two implementation 
lambdahat - lambdahat_R 

par(mfrow=c(1,1))
curve(exp(sapply(x, Qfun, G = G)), from = lambdahat-1, to = lambdahat+1)
abline(v = lambdahat, col='red')
abline(v = lambdahat_R, col='blue')

# multi-dimensional problem
N <- 10 # number of observations
m <- 3 # number of dimensions
y <- rnorm(N*m)
X <- matrix(rep(1,N))
G <- matrix(rnorm(N*m), m, N)

# optimization in C++
lambda0 <- rnorm(m)
lambdahat <- lambdaNR(y = y, X = X, G = t(G), lambda0 = lambda0)

# by implementation in R
lambdahat_Rout <- lambdaNR_R(G = t(G), lambda0 = lambda0)
lambdahat_R <- t(lambdahat_Rout$lambda) # output is a row vector

# difference of the two implementation 
lambdahat - lambdahat_R 

# derivative test 
require(numDeriv)

# grad Q(lambdahat) = 0 (sometimes unreliable)
Qf <- function(lambda) Qfun(lambda, G)
grad(Qf, lambdahat)
grad(Qf, lambdahat_R)

# visual mode check (more reliable)
mle.check(loglik = Qf, theta.mle = lambdahat)

# visual mode check (more reliable)
mle.check(loglik = Qf, theta.mle = lambdahat_R)
