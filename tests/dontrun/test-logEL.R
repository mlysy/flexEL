#--- check that logELmean.R is working properly ---------------------------------

require(bayesEL)
source("../bayesEL/tests/el-utils.R")

# Plot the EL likelihood EL(mu) for the single constraint E[Y - mu] = 0

# dimensions
n <- 100 # number of observations
numpoints <- 100

#---- mean reg: p = 1 (only intercept) ----
m <- 2 # number of estimating equations [always 2 now, mean = 0 and var = 1]
p <- 1 # dimension of beta
mu0 <- rnorm(p) # true parameter value
mu0

# y <- rnorm(n, mean = mu0) # data
# y <- rchisq(n, df=3)/sqrt(6)
# mu0 <- 3/sqrt(6)
# mean(y)

# gird plot
mu.seq <- seq(-.5+mu0,.5+mu0,length.out = numpoints)
logel.seq <- rep(NA,numpoints)
X <- matrix(rep(1,p*n),p,n) # each col of X is one observation
y <- X * mu0 + rnorm(n) # with N(0,1) error term
mean(y)
lambda0 <- rnorm(m)
for (ii in 1:numpoints) {
      logel.seq[ii] <- logELmean(n, m, y, X, mu.seq[ii],
                                 lambda0, max_iter = 100, eps = 1e-7)
}
logelmode <- plotEL(mu.seq, logel.seq, mu0, mean(y), expression(mu))

# sample from posterior
nsamples <- 10000
nburn <- 2000
beta0 <- rnorm(p) # current prior of beta N(0,1)
sigs <- rep(0.1,m)
system.time(
  mu_chain <- MeanReg_post(n, m, y, X, lambda0, nsamples, nburn, beta0, sigs)
)
plot(mu_chain[1,], xlab = 'mu', ylab = 'EL', type='l')

# overlay gird plot to histogram
hist(mu_chain[1,],breaks=50,freq=FALSE,
     xlab = expression(mu),main='')
lines(mu.seq,norm_pdf(logel.seq,mu.seq),
      cex=0.1, col = 'red', type='l')
abline(v=logelmode, col='red')
abline(v=mean(mu_chain), col='blue')
legend('topright',legend=c(expression('grid plot & mode'),
                           expression('sample mean')),
       lty = c(1,1), col = c('red','blue'), cex = 0.6)

#---- mean reg: p = 2 (one intercept and one covariate) ----
m <- 2
p <- 2
X <- rbind(rep(1,n),rnorm(n))
beta0 <- rnorm(p)
y <- t(X) %*% beta0 + rnorm(n) # with N(0,1) error term
mu.seq.b1 <- seq(-.5+beta0[1],.5+beta0[1],length.out = numpoints) # TODO: grid points range matters ?
mu.seq.b2 <- seq(-.5+beta0[2],.5+beta0[2],length.out = numpoints)
mu.seq <- rbind(mu.seq.b1, mu.seq.b2)
lambda0 <- rnorm(m)
logel.seq <- matrix(rep(NA,2*numpoints),2,numpoints)
for (ii in 1:numpoints) {
      logel.seq[1,ii] <- logELmean(n, m, y, X, c(mu.seq[1,ii], beta0[2]), lambda0)
      logel.seq[2,ii] <- logELmean(n, m, y, X, c(beta0[1], mu.seq[2,ii]), lambda0)
}
logelmode1 <- plotEL(mu.seq[1,], logel.seq[1,], beta0[1], NA, expression(beta[0]))
logelmode2 <- plotEL(mu.seq[2,], logel.seq[2,], beta0[2], NA, expression(beta[1]))

# sample from posterior
nsamples <- 10000
nburn <- 2000
sigs <- rep(0.1,m)
betaInit <- rnorm(p)
system.time(
  beta_chain <- MeanReg_post(n, m, y, X, lambda0, nsamples, nburn, betaInit, sigs)
)
plot(beta_chain[1,], xlab = 'beta1', ylab = 'EL', type='l')
plot(beta_chain[2,], xlab = 'beta2', ylab = 'EL', type='l')

# overlay gird plot with histogram 
hist(beta_chain[1,],breaks=50,freq=FALSE,
     xlab = expression(beta[0]),main='')
lines(mu.seq.b1,norm_pdf(logel.seq[1,],mu.seq.b1),
      cex=0.1, col = 'red', type='l')
abline(v=logelmode1, col='red')
abline(v=mean(beta_chain[1,]), col='blue')
legend('topright',legend=c(expression('grid plot & mode'),
                           expression('sample mean')),
       lty = c(1,1), col = c('red','blue'), cex = 0.6)

hist(beta_chain[2,],breaks=50,freq=FALSE,
     xlab = expression(beta[1]),main='')
lines(mu.seq.b2,norm_pdf(logel.seq[2,],mu.seq.b2),
      cex=0.1, col = 'red', type='l')
abline(v=logelmode2, col='red')
abline(v=mean(beta_chain[2,]), col='blue')
legend('topright',legend=c(expression('grid plot & mode'),
                           expression('sample mean')),
       lty = c(1,1), col = c('red','blue'), cex = 0.6)

#---- quantile reg:  p = 1 (only intercept) ----
m <- 1 # for quant reg m = p
mu0 <- rnorm(m) # true parameter value
mu0

# gird plot
mu.seq <- seq(-.5+mu0,.5+mu0,length.out = numpoints) # TODO: grid points range matters ?
logel.seq <- rep(NA,numpoints)
X <- matrix(rep(1,m*n),m,n) # each col of X is one observation
y <- X * mu0 + rnorm(n) # with N(0,1) error term
mean(y)
lambda0 <- rnorm(m)
alpha <- 0.5
logel.seq <- rep(NA,numpoints)
for (ii in 1:numpoints) {
  logel.seq[ii] <- logELquant(n, m, X, y, mu.seq[ii], alpha, lambda0)
}
logelmode <- plotEL(mu.seq, logel.seq, mu0, mean(y), expression(mu))

# sample from posteriro
nsamples <- 10000
nburn <- 2000
sigs <- rep(0.1,m)
betaInit <- rnorm(m)
system.time(
  beta_chain <- QuantReg_post(n, m, y, X, alpha, lambda0, nsamples, nburn, betaInit, sigs)
)
plot(beta_chain[1,], xlab = expression(mu), ylab = 'logEL', type='l')

# overlay gird plot to histogram
hist(beta_chain[1,],breaks=50,freq=FALSE,
     xlab = expression(mu),main='')
lines(mu.seq,norm_pdf(logel.seq,mu.seq),
      cex=0.1, col = 'red', type='l')
abline(v=logelmode, col='red')
abline(v=mean(beta_chain[1,]), col='blue')
legend('topright',legend=c(expression('grid plot & mode'),
                           expression('sample mean')),
       lty = c(1,1), col = c('red','blue'), cex = 0.6)

#---- quantile reg:  p = 2 (one intercept and one covariate) ----

# Plot the EL for quantile regression with constraint E[phi(Y - mu)] = 0
m <- 2 #  for quant reg m = p
X <- rbind(rep(1,n),rnorm(n))
beta0 <- rnorm(m)
beta0
y <- t(X) %*% beta0 + rnorm(n)
mu.seq.b1 <- seq(beta0[1]-.5,beta0[1]+.5,length.out = numpoints)
mu.seq.b2 <- seq(beta0[2]-.5,beta0[2]+.5,length.out = numpoints)
mu.seq <- rbind(mu.seq.b1,mu.seq.b2)
logel.seq <- matrix(rep(NA,2*numpoints),2,numpoints)
lambda0 <- rnorm(m)
alpha <- 0.5 # Seems extreme quantile e.g. 0.9 is much harder
for (ii in 1:numpoints) {
      logel.seq[1,ii] <- logELquant(n, m, X, y, c(mu.seq[1,ii],beta0[2]), alpha, lambda0)
      logel.seq[2,ii] <- logELquant(n, m, X, y, c(beta0[1],mu.seq[2,ii]), alpha, lambda0)
}
logelmode1 <- plotEL(mu.seq[1,], logel.seq[1,], beta0[1], NA, expression(beta[0]))
logelmode2 <- plotEL(mu.seq[2,], logel.seq[2,], beta0[2], NA, expression(beta[1]))

# sample from posterior
nsamples <- 10000
nburn <- 2000
sigs <- rep(0.1,m)
betaInit <- rnorm(m)
system.time(
  beta_chain <- QuantReg_post(n, m, y, X, alpha, lambda0, nsamples, nburn, betaInit, sigs)
)
plot(beta_chain[1,], xlab = 'beta1', ylab = 'EL', type='l')
plot(beta_chain[2,], xlab = 'beta2', ylab = 'EL', type='l')

# overlay gird plot to histogram
hist(beta_chain[1,],breaks=50,freq=FALSE,
     xlab = expression(beta[0]),main='')
lines(mu.seq.b1,norm_pdf(logel.seq[1,],mu.seq.b1),
      cex=0.1, col = 'red', type='l')
abline(v=logelmode1, col='red')
abline(v=mean(beta_chain[1,]), col='blue')
legend('topright',legend=c(expression('grid plot & mode'),
                           expression('sample mean')),
       lty = c(1,1), col = c('red','blue'), cex = 0.6)

hist(beta_chain[2,],breaks=50,freq=FALSE,
     xlab = expression(beta[1]),main='')
lines(mu.seq.b2,norm_pdf(logel.seq[2,],mu.seq.b2),
      cex=0.1, col = 'red', type='l')
abline(v=logelmode2, col='red')
abline(v=mean(beta_chain[2,]), col='blue')
legend('topright',legend=c(expression('grid plot & mode'),
                           expression('sample mean')),
       lty = c(1,1), col = c('red','blue'), cex = 0.6)

