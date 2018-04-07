library(bayesEL)
source("../testthat/el-utils.R")
# tests for quantile regression 
# ---- 1-d problem ----
n <- 500
p <- 1
mu <- 1
alpha <- 0.75
X <- matrix(rep(1,p*n), n, p) # each row of X is one observation
eps <- rnorm(n) # N(0,1) error term
quantile(eps, alpha)
y <- c(X * mu) + eps
mu0 <- mu + qnorm(alpha, lower.tail = TRUE) # true param assuming alpha-tail of eps is 0
mu0

# gird plot
numpoints <- 100
mu.seq <- seq(-.5+mu0,.5+mu0,length.out = numpoints)
logel.seq <- rep(NA,numpoints)
for (ii in 1:numpoints) {
  G <- qr.evalG(y,X,alpha,mu.seq[ii])
  logel.seq[ii] <- logEL(G = G)
}
logelmode <- plotEL(mu.seq, logel.seq, mu0, quantile(y,alpha), expression(mu))
logelmode

nsamples <- 20000
nburn <- 5000
# betaInit <- mu0
betaInit <- rnorm(length(mu0), mean = mu0, sd = 1) # TODO: if init far away, problematic ..
betaInit
sigs <- rep(0.01,1)
mu_chain <- qr.post(y, X, alpha, nsamples, nburn, betaInit, sigs)
plot(mu_chain[1,], xlab = 'mu', ylab = 'EL', type='l')
# overlay gird plot to histogram
hist(mu_chain[1,],breaks=50,freq=FALSE,
     xlab = expression(mu), main='')
lines(mu.seq, norm_pdf(logel.seq, mu.seq),
      cex=0.1, col = 'red', type='l')
abline(v=logelmode, col='red')
abline(v=mean(mu_chain), col='blue')
legend('topright',legend=c(expression('grid plot & mode'),
                           expression('sample mean')),
       lty = c(1,1), col = c('red','blue'), cex = 0.6)

# ---- 2-d problem (1 intercept, 1 slope) ----
n <- 500
p <- 2
alpha <- 0.75
# X0 <- matrix(rep(1,n),n,1)
X1 <- matrix(seq(-2,2,length.out = n),n,1)
X1 <- matrix(rnorm(n),n,1)
X <- cbind(X0,X1)
eps <- rnorm(n) # N(0,1) error term
beta_intercept <- qnorm(alpha, lower.tail = TRUE)
beta_slope <- 1.5
y <- c(X1 %*% beta_slope) + eps 
plot(X1,y,cex=0.3)
beta0 <- c(beta_intercept, beta_slope)

# grid plot
numpoints <- 100
beta1.seq <- seq(beta0[1]-.5,beta0[1]+.5,length.out = numpoints)
beta2.seq <- seq(beta0[2]-.5,beta0[2]+.5,length.out = numpoints)
logel.seq <- matrix(rep(NA,2*numpoints),2,numpoints)
for (ii in 1:numpoints) {
  G <- qr.evalG(y,X,alpha,c(beta1.seq[ii],beta0[2]))
  logel.seq[1,ii] <- logEL(G)
  G <- qr.evalG(y,X,alpha,c(beta0[1],beta2.seq[ii]))
  logel.seq[2,ii] <- logEL(G)
}
logelmode1 <- plotEL(beta1.seq, logel.seq[1,], beta0[1], NA, expression(beta[0]))
logelmode2 <- plotEL(beta2.seq, logel.seq[2,], beta0[2], NA, expression(beta[1]))

nsamples <- 20000
nburn <- 5000
betaInit <- beta0
sigs <- rep(0.2,2)
system.time(
  beta_chain <- qr.post(y, X, alpha, nsamples, nburn, betaInit, sigs)
)
plot(beta_chain[1,], xlab = expression(beta[0]), ylab = 'EL', type='l')
plot(beta_chain[2,], xlab = expression(beta[1]), ylab = 'EL', type='l')
# overlay gird plot to histogram
# intercept
hist(beta_chain[1,],breaks=50,freq=FALSE,
     xlab = expression(beta[0]),main='')
lines(beta1.seq, norm_pdf(logel.seq[1,], beta1.seq),
      cex=0.1, col = 'red', type='l')
abline(v=logelmode1, col='red')
abline(v=mean(beta_chain[1,]), col='blue')
legend('topright',legend=c(expression('grid plot & mode'),
                           expression('sample mean')),
       lty = c(1,1), col = c('red','blue'), cex = 0.6)
# slope
hist(beta_chain[2,],breaks=50,freq=FALSE,
     xlab = expression(beta[1]),main='')
lines(beta2.seq, norm_pdf(logel.seq[2,], beta2.seq),
      cex=0.1, col = 'red', type='l')
abline(v=logelmode2, col='red')
abline(v=mean(beta_chain[2,]), col='blue')
legend('topright',legend=c(expression('grid plot & mode'),
                           expression('sample mean')),
       lty = c(1,1), col = c('red','blue'), cex = 0.6)

# ---- 3-d problem (1 intercept, 2 slope) ---- 

