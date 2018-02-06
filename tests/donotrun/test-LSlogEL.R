# check that MeanRegLS_logEL.R is working properly
require(bayesEL)
source("../bayesEL/tests/testthat/el-utils.R")

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

G_R <- LSevalG(y, X, theta0)
