# tests for quantile regression postSampler with adapt MCMC
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
  # G <- qr.evalG(y,X,alpha,mu.seq[ii])
  G <- qr.evalG.smooth_R(y,X,alpha,mu.seq[ii])
  omegas <- omega.hat(G)
  logel.seq[ii] <- logEL(omegas)
}
logelmode <- plotEL(mu.seq, logel.seq, mu0, quantile(y,alpha), expression(mu))
logelmode

# # plot two figures side by side (compare the smoothing effect)
# par(mfrow=c(1,2))


nsamples <- 20000
nburn <- 5000
library(quantreg)
betaInit <- c(rq(y ~ 1, tau = alpha, method = 'fn')$coefficients)
betaInit
sigs <- rep(0.25,1)
rvDoMcmc <- c(1)
qrout <- qr.post_adapt(y, X, alpha, nsamples, nburn, betaInit, sigs, rvDoMcmc)
mu_chain <- qrout$beta_chain
mu_paccept <- qrout$paccept 
mu_paccept
plot(mu_chain[1,], xlab = 'mu', ylab = 'EL', type='l')
# overlay gird plot to histogram
hist(mu_chain[1,],breaks=50,freq=FALSE,
     xlab = expression(mu), main='')
lines(mu.seq, norm_pdf(logel.seq, mu.seq),
      cex=0.1, col = 'red', type='l')
abline(v=mu0, col='red')
abline(v=mean(mu_chain), col='blue')
legend('topright',legend=c(expression('true param'),
                           expression('sample mean')),
       lty = c(1,1), col = c('red','blue'), cex = 0.6)

# ---- 2-d problem (1 intercept, 1 slope) ----
n <- 200
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
  G <- qr.evalG(y,X,alpha,c(beta1.seq[ii],beta0[2]))
  omegas <- omega.hat(G)
  logel.seq[1,ii] <- logEL(omegas)
  G <- qr.evalG(y,X,alpha,c(beta0[1],beta2.seq[ii]))
  omegas <- omega.hat(G)
  logel.seq[2,ii] <- logEL(omegas)
}
par(mfrow=c(1,2))
logelmode1 <- plotEL(beta1.seq, logel.seq[1,], beta0[1], NA, expression(beta[0]),
                     legend.loc = 'topright', cex.lab = 1.5, cex.axis = 1.2, cex=1.2)
logelmode2 <- plotEL(beta2.seq, logel.seq[2,], beta0[2], NA, expression(beta[1]),
                     legend.loc = 'topleft', cex.lab = 1.5, cex.axis = 1.2, cex=1.2)

# calculate marginal posterior
Beta.seq <- as.matrix(expand.grid(beta1.seq, beta2.seq))
logel.mat <- apply(Beta.seq, 1, function(bb) {
  G <- qr.evalG(y,X,alpha,c(bb[1],bb[2]))
  logEL(G)
})
logel.mat <- matrix(logel.mat, numpoints, numpoints)
el.mat <- exp(logel.mat - max(logel.mat))
logel.marg <- log(cbind(beta1 = rowSums(el.mat), beta2 = colSums(el.mat)))

library(quantreg)
betaInit <- c(rq(y ~ X1, tau = alpha, method = 'fn')$coefficients)
betaInit
nsamples <- 20000
nburn <- 3000
sigs <- rep(0.2,2)
rvDoMcmc <- c(0,1)
system.time(
  qrout <- qr.post_adapt(y, X, alpha, nsamples, nburn, c(beta0[1],betaInit[2]), sigs, rvDoMcmc)
)
beta_chain <- qrout$beta_chain
beta_paccept <- qrout$paccept
beta_paccept
plot(beta_chain[1,], xlab = expression(beta[0]), ylab = 'EL', type='l')
plot(beta_chain[2,], xlab = expression(beta[1]), ylab = 'EL', type='l')
# overlay gird plot to histogram
# intercept
hist(beta_chain[1,],breaks=50,freq=FALSE,
     xlab = expression(beta[0]),main='')
# marginal line
lines(beta1.seq, norm_pdf(logel.marg[,1], beta1.seq),
      cex=0.1, col = 'red', type='l')
# conditional line
lines(beta1.seq, norm_pdf(logel.seq[1,], beta1.seq),
      cex=0.1, col = 'red', type='l')
abline(v=beta0[1], col='red')
abline(v=mean(beta_chain[1,]), col='blue')
legend('topright',legend=c(expression('true param'),
                           expression('sample mean')),
       lty = c(1,1), col = c('red','blue'), cex = 0.6)
# slope
hist(beta_chain[2,],breaks=50,freq=FALSE,
     xlab = expression(beta[1]),main='')
# marginal line
lines(beta2.seq, norm_pdf(logel.marg[,2], beta2.seq),
      cex=0.1, col = 'red', type='l')
# conditional line
lines(beta2.seq, norm_pdf(logel.seq[2,], beta2.seq),
      cex=0.1, col = 'red', type='l')
abline(v=beta0[2], col='red')
abline(v=mean(beta_chain[2,]), col='blue')
legend('topright',legend=c(expression('true param'),
                           expression('sample mean')),
       lty = c(1,1), col = c('red','blue'), cex = 0.6)
