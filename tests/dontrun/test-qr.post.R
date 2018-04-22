# tests for quantile regression postSampler 
library(bayesEL)
source("../testthat/el-utils.R")
# ---- 1-d problem ----
n <- 500
mu <- 1
alpha <- 0.75
X <- matrix(rep(1,n), n, 1) # each row of X is one observation
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
library(quantreg)
betaInit <- c(rq(y ~ 1, tau = alpha, method = 'fn')$coefficients)
betaInit
sigs <- rep(0.25,1)
qrout <- qr.post(y, X, alpha, nsamples, nburn, betaInit, sigs)
mu_chain <- qrout$Beta_chain
mu_paccept <- qrout$paccept 
mu_paccept
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

# ---- 1-d problem with 2 quantile levels ---- 
n <- 100
mu <- 1
alphas <- c(0.5,0.75)
X <- matrix(rep(1,n), n, 1) # each row of X is one observation
eps <- rnorm(n) # N(0,1) error term
quantile(eps, alphas[1])
quantile(eps, alphas[2])
y <- c(X * mu) + eps
mu01 <- mu + qnorm(alphas[1], lower.tail = TRUE) # true param assuming alpha-tail of eps is 0
mu01
mu02 <- mu + qnorm(alphas[2], lower.tail = TRUE) # true param assuming alpha-tail of eps is 0
mu02

# gird plot
numpoints <- 100
mu.seq1 <- seq(-.5+mu01,.5+mu01,length.out = numpoints)
logel.seq1 <- rep(NA,numpoints)
for (ii in 1:numpoints) {
  G <- qr.evalG(y,X,alphas[1],mu.seq1[ii])
  logel.seq1[ii] <- logEL(G = G)
}
logelmode1 <- plotEL(mu.seq1, logel.seq1, mu01, quantile(y,alphas[1]), expression(mu))
logelmode1

numpoints <- 100
mu.seq2 <- seq(-.5+mu02,.5+mu02,length.out = numpoints)
logel.seq2 <- rep(NA,numpoints)
for (ii in 1:numpoints) {
  G <- qr.evalG(y,X,alphas[2],mu.seq2[ii])
  logel.seq2[ii] <- logEL(G = G)
}
logelmode2 <- plotEL(mu.seq2, logel.seq2, mu02, quantile(y,alphas[2]), expression(mu))
logelmode2

nsamples <- 20000
nburn <- 3000
library(quantreg)
betaInit1 <- c(rq(y ~ 1, tau = alphas[1], method = 'fn')$coefficients)
betaInit1
betaInit2 <- c(rq(y ~ 1, tau = alphas[2], method = 'fn')$coefficients)
betaInit2
BetaInit <- cbind(betaInit1,betaInit2)
# BetaInit <- cbind(mu01,mu02)
sigs1 <- rep(0.22,1)
sigs2 <- rep(0.15,1)
Sigs <- cbind(sigs1,sigs2)
qrout <- qr.post(y, X, alphas, nsamples, nburn, BetaInit, Sigs)

mu_chain <- qrout$Beta_chain
mu_paccept <- qrout$paccept 
mu_paccept
plot(mu_chain[1,], xlab = 'mu', ylab = 'EL', type='l')
plot(mu_chain[2,], xlab = 'mu', ylab = 'EL', type='l')

# overlay gird plot to histogram
hist(mu_chain[1,],breaks=50,freq=FALSE,
     xlab = expression(mu), main='')
lines(mu.seq1, norm_pdf(logel.seq1, mu.seq1),
      cex=0.1, col = 'red', type='l')
abline(v=logelmode1, col='red')
abline(v=mean(mu_chain[1,]), col='blue')
legend('topright',legend=c(expression('grid plot & mode'),
                           expression('sample mean')),
       lty = c(1,1), col = c('red','blue'), cex = 0.6)

hist(mu_chain[2,],breaks=50,freq=FALSE,
     xlab = expression(mu), main='')
lines(mu.seq2, norm_pdf(logel.seq2, mu.seq2),
      cex=0.1, col = 'red', type='l')
abline(v=logelmode2, col='red')
abline(v=mean(mu_chain[2,]), col='blue')
legend('topright',legend=c(expression('grid plot & mode'),
                           expression('sample mean')),
       lty = c(1,1), col = c('red','blue'), cex = 0.6)

# ---- 2-d problem (1 intercept, 1 slope) ----
n <- 500
p <- 2
alpha <- 0.75
X0 <- matrix(rep(1,n),n,1)
# X1 <- matrix(seq(-2,2,length.out = n),n,1)
X1 <- matrix(rnorm(n),n,1)
X <- cbind(X0,X1)
eps <- rnorm(n) # N(0,1) error term
beta_intercept <- qnorm(alpha, lower.tail = TRUE)
beta_slope <- 1.5
y <- c(X1 %*% beta_slope) + eps 
plot(X1,y,cex=0.3)
beta0 <- c(beta_intercept, beta_slope)
beta0

# grid plot (conditionals: beta1|beta2 and beta2|beta1)
numpoints <- 100
beta1.seq <- seq(beta0[1]-.5,beta0[1]+.5,length.out = numpoints)
beta2.seq <- seq(beta0[2]-.5,beta0[2]+.5,length.out = numpoints)
logel.seq <- matrix(rep(NA,2*numpoints),2,numpoints)
for (ii in 1:numpoints) {
  # G <- qr.evalG(y,X,alpha,c(beta1.seq[ii],beta0[2]))
  G <- qr.evalG(y,X,alpha,c(beta1.seq[ii],beta0[2]))
  logel.seq[1,ii] <- logEL(G)
  G <- qr.evalG(y,X,alpha,c(beta0[1],beta2.seq[ii]))
  logel.seq[2,ii] <- logEL(G)
}
logelmode1 <- plotEL(beta1.seq, logel.seq[1,], beta0[1], NA, expression(beta[0]))
logelmode2 <- plotEL(beta2.seq, logel.seq[2,], beta0[2], NA, expression(beta[1]))

# calculate marginal posterior
Beta.seq <- as.matrix(expand.grid(beta1.seq, beta2.seq))
logel.mat <- apply(Beta.seq, 1, function(bb) {
  G <- qr.evalG(y,X,alpha,matrix(c(bb[1],bb[2]),2,1))
  logEL(G)
})
logel.mat <- matrix(logel.mat, numpoints, numpoints)
el.mat <- exp(logel.mat - max(logel.mat))
logel.marg <- log(cbind(beta1 = rowSums(el.mat), beta2 = colSums(el.mat)))

nsamples <- 20000
nburn <- 5000
library(quantreg)
betaInit <- c(rq(y ~ X1, tau = alpha, method = 'fn')$coefficients)
betaInit
sigs <- rep(0.2,2)
system.time(
  qrout <- qr.post(y, X, alpha, nsamples, nburn, betaInit, sigs)
)
beta_chain <- qrout$Beta_chain
beta_paccept <- qrout$paccept
beta_paccept
plot(beta_chain[1,], xlab = expression(beta[0]), ylab = 'EL', type='l')
plot(beta_chain[2,], xlab = expression(beta[1]), ylab = 'EL', type='l')
# overlay gird plot to histogram
# intercept
hist(beta_chain[1,],breaks=50,freq=FALSE,
     xlab = expression(beta[0]),main='')
lines(beta1.seq, norm_pdf(logel.marg[,1], beta1.seq),
      cex=0.1, col = 'red', type='l')
abline(v=mean(beta_chain[1,]), col='blue')
legend('topright',legend=c(expression('grid plot'),
                           expression('sample mean')),
       lty = c(1,1), col = c('red','blue'), cex = 0.6)
# slope
hist(beta_chain[2,],breaks=50,freq=FALSE,
     xlab = expression(beta[1]),main='')
lines(beta2.seq, norm_pdf(logel.marg[,2], beta2.seq),
      cex=0.1, col = 'red', type='l')
abline(v=mean(beta_chain[2,]), col='blue')
legend('topright',legend=c(expression('grid plot & mode'),
                           expression('sample mean')),
       lty = c(1,1), col = c('red','blue'), cex = 0.6)

# ---- 2-d problem with 2 quantile levels (1 intercept, 1 slope) ----
n <- 500
alphas <- c(0.75,0.9)
X0 <- matrix(rep(1,n),n,1)
# X1 <- matrix(seq(-2,2,length.out = n),n,1)
X1 <- matrix(rnorm(n),n,1)
X <- cbind(X0,X1)
eps <- rnorm(n) # N(0,1) error term

beta_slope <- 1.5
beta_intercept1 <- qnorm(alphas[1], lower.tail = TRUE)
beta01 <- c(beta_intercept1, beta_slope)
beta01
beta_intercept2 <- qnorm(alphas[2], lower.tail = TRUE)
beta02 <- c(beta_intercept2, beta_slope)
beta02
y <- c(X1 %*% beta_slope) + eps 
plot(X1,y,cex=0.3)

# grid plot (conditionals: beta1|beta2 and beta2|beta1)
numpoints <- 100
beta1.seq1 <- seq(beta01[1]-.5,beta01[1]+.5,length.out = numpoints)
beta2.seq1 <- seq(beta01[2]-.5,beta01[2]+.5,length.out = numpoints)
logel.seq1 <- matrix(rep(NA,2*numpoints),2,numpoints)
for (ii in 1:numpoints) {
  G <- qr.evalG(y,X,alphas[1],c(beta1.seq1[ii],beta01[2]))
  logel.seq1[1,ii] <- logEL(G)
  G <- qr.evalG(y,X,alphas[1],c(beta01[1],beta2.seq1[ii]))
  logel.seq1[2,ii] <- logEL(G)
}
logelmode1 <- plotEL(beta1.seq1, logel.seq1[1,], beta01[1], NA, expression(beta[0]))
logelmode2 <- plotEL(beta2.seq1, logel.seq1[2,], beta01[2], NA, expression(beta[1]))

# calculate marginal posterior
Beta.seq1 <- as.matrix(expand.grid(beta1.seq1, beta2.seq1))
logel.mat1 <- apply(Beta.seq1, 1, function(bb) {
  G <- qr.evalG(y,X,alphas[1],matrix(c(bb[1],bb[2]),2,1))
  logEL(G)
})
logel.mat1 <- matrix(logel.mat1, numpoints, numpoints)
el.mat1 <- exp(logel.mat1 - max(logel.mat1))
logel.marg1 <- log(cbind(beta1 = rowSums(el.mat1), beta2 = colSums(el.mat1)))

# grid plot (conditionals: beta1|beta2 and beta2|beta1)
numpoints <- 100
beta1.seq2 <- seq(beta02[1]-.5,beta02[1]+.5,length.out = numpoints)
beta2.seq2 <- seq(beta02[2]-.5,beta02[2]+.5,length.out = numpoints)
logel.seq2 <- matrix(rep(NA,2*numpoints),2,numpoints)
for (ii in 1:numpoints) {
  G <- qr.evalG(y,X,alphas[2],c(beta1.seq2[ii],beta02[2]))
  logel.seq2[1,ii] <- logEL(G)
  G <- qr.evalG(y,X,alphas[2],c(beta02[1],beta2.seq2[ii]))
  logel.seq2[2,ii] <- logEL(G)
}
logelmode3 <- plotEL(beta1.seq2, logel.seq2[1,], beta02[1], NA, expression(beta[0]))
logelmode4 <- plotEL(beta2.seq2, logel.seq2[2,], beta02[2], NA, expression(beta[1]))

# calculate marginal posterior
Beta.seq2 <- as.matrix(expand.grid(beta1.seq2, beta2.seq2))
logel.mat2 <- apply(Beta.seq2, 1, function(bb) {
  G <- qr.evalG(y,X,alphas[2],matrix(c(bb[1],bb[2]),2,1))
  logEL(G)
})
logel.mat2 <- matrix(logel.mat2, numpoints, numpoints)
el.mat2 <- exp(logel.mat2 - max(logel.mat2))
logel.marg2 <- log(cbind(beta1 = rowSums(el.mat2), beta2 = colSums(el.mat2)))

nsamples <- 20000
nburn <- 5000
library(quantreg)
betaInit1 <- c(rq(y ~ X1, tau = alphas[1], method = 'fn')$coefficients)
betaInit1
betaInit2 <- c(rq(y ~ X1, tau = alphas[2], method = 'fn')$coefficients)
betaInit2
BetaInit <- cbind(betaInit1,betaInit2)
# Sigs <- cbind(c(0.27,0.23),c(0.27,0.26))
Sigs <- cbind(c(0.13,0.12),c(0.18,0.18))
system.time(
  qrout <- qr.post(y, X, alphas, nsamples, nburn, BetaInit, Sigs)
)
beta_chain <- qrout$Beta_chain
beta_paccept <- qrout$paccept
beta_paccept
plot(beta_chain[1,], xlab = expression(beta[0]), ylab = 'EL', type='l')
plot(beta_chain[2,], xlab = expression(beta[1]), ylab = 'EL', type='l')
plot(beta_chain[3,], xlab = expression(beta[0]), ylab = 'EL', type='l')
plot(beta_chain[4,], xlab = expression(beta[1]), ylab = 'EL', type='l')
# overlay gird plot to histogram
# intercept
hist(beta_chain[1,],breaks=50,freq=FALSE,
     xlab = expression(beta[0]),main='')
lines(beta1.seq1, norm_pdf(logel.marg1[,1], beta1.seq1),
      cex=0.1, col = 'red', type='l')
abline(v=beta01[1], col='red')
abline(v=mean(beta_chain[1,]), col='blue')
legend('topright',legend=c(expression('grid plot & true param'),
                           expression('sample mean')),
       lty = c(1,1), col = c('red','blue'), cex = 0.6)
# slope
hist(beta_chain[2,],breaks=50,freq=FALSE,
     xlab = expression(beta[1]),main='')
lines(beta2.seq1, norm_pdf(logel.marg1[,2], beta2.seq1),
      cex=0.1, col = 'red', type='l')
abline(v=beta01[2], col='red')
abline(v=mean(beta_chain[2,]), col='blue')
legend('topright',legend=c(expression('grid plot & true param'),
                           expression('sample mean')),
       lty = c(1,1), col = c('red','blue'), cex = 0.6)
# intercept
hist(beta_chain[3,],breaks=50,freq=FALSE,
     xlab = expression(beta[0]),main='')
lines(beta1.seq2, norm_pdf(logel.marg2[,1], beta1.seq2),
      cex=0.1, col = 'red', type='l')
abline(v=beta02[1], col='red')
abline(v=mean(beta_chain[3,]), col='blue')
legend('topright',legend=c(expression('grid plot & true param'),
                           expression('sample mean')),
       lty = c(1,1), col = c('red','blue'), cex = 0.6)
# slope
hist(beta_chain[4,],breaks=50,freq=FALSE,
     xlab = expression(beta[1]),main='')
lines(beta2.seq2, norm_pdf(logel.marg2[,2], beta2.seq2),
      cex=0.1, col = 'red', type='l')
abline(v=beta02[2], col='red')
abline(v=mean(beta_chain[4,]), col='blue')
legend('topright',legend=c(expression('grid plot & true param'),
                           expression('sample mean')),
       lty = c(1,1), col = c('red','blue'), cex = 0.6)

