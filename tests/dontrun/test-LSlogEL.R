# check that MeanRegLS_logEL.R is working properly
require(bayesEL)
#source("../bayesEL/tests/testthat/el-utils.R")
source("../testthat/el-utils.R")

# dimensions
n <- 5 # number of observations
numpoints <- 100

#---- mean reg: p = 1 (only intercept) ----
p <- 1
X <- matrix(rep(1,n),1,n)
# beta0 <- rnorm(p)
# gamma0 <- rnorm(p)
beta0 <- 2
gamma0 <- -0.01
theta0 <- c(beta0,gamma0)
y <- t(X) %*% beta0 + exp(t(X) %*% gamma0)*rnorm(n)
plot(y)
beta.seq <- seq(-.5+beta0, .5+beta0, length.out = numpoints)
gamma.seq <- seq(-.5+gamma0, .5+gamma0, length.out = numpoints)
theta.seq <- rbind(beta.seq, gamma.seq)
logel.seq <- matrix(rep(NA,2*numpoints),2,numpoints)
for (ii in 1:numpoints) {
    logel.seq[1,ii] <- mrls.logel(y, X, c(theta.seq[1,ii], gamma0))
    logel.seq[2,ii] <- mrls.logel(y, X, c(beta0, theta.seq[2,ii]))
}
logelmode1 <- plotEL(beta.seq, logel.seq[1,], beta0, NA, expression(beta))
logelmode2 <- plotEL(gamma.seq, logel.seq[2,], gamma0, NA, expression(gamma))

#--- test G function -----------------------------------------------------------

nreps <- 100
max.abs <- replicate(nreps, expr = {
  nObs <- sample(10:100, 1)
  nBeta <- sample(2:10, 1)
  nGamma <- sample(2:10, 1)
  X <- rMNorm(nObs, nBeta)
  beta <- rnorm(nBeta)
  Z <- rMNorm(nObs, nGamma)
  gamma <- rnorm(nGamma)
  y <- X %*% beta + exp(Z %*% gamma) * rnorm(nObs)
  y <- c(y)
  G_R <- mrls.evalG_R(y, X, Z, beta, gamma)
  G_cpp <- mrls.G(y, X, Z, beta, gamma)
  abs(G_R - t(G_cpp))
})
max(sapply(max.abs, max))

# test loglikelihood

nObs <- sample(10:100, 1)
nBeta <- sample(2:10, 1)
nGamma <- sample(2:10, 1)
X <- rMNorm(nObs, nBeta)
beta0 <- rnorm(nBeta)
Z <- rMNorm(nObs, nGamma)
gamma0 <- rnorm(nGamma)
y <- X %*% beta0 + exp(Z %*% gamma0) * rnorm(nObs)
y <- c(y)

mrls.logel(y, X, Z, beta0, gamma0)
mrls.logel_R(y, X, Z, beta0, gamma0)


nreps <- 100
Theta <- rMNorm(nreps, nBeta+nGamma)

apply(Theta, 1, function(theta) {
  beta <- theta[1:nBeta]
  gamma <- theta[nBeta + 1:nGamma]
  c(mrls.logel(y, X, Z, beta, gamma), mrls.logel_R(y, X, Z, beta, gamma))
})

