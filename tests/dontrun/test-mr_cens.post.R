# tests for mean regression postSampler 
library(bayesEL)
source("../testthat/el-utils.R")
source("../testthat/el-rfuns.R")
source("../testthat/el-model.R")

# ---- 1-d problem ----
n <- 30
mu0 <- rnorm(1,0,1)
X <- matrix(rnorm(n,0,1),n,1) # each row of X is one observation
# X <- matrix(rep(1,n), n, 1)
eps <- rnorm(n) # N(0,1) error term
yy <- c(X * mu0) + eps

# random censoring
cc <- rnorm(n,mean=1.5*abs(yy),sd=1)
deltas <- yy<=cc
y <- yy
sum(1-deltas)/n
y[as.logical(1-deltas)] <- cc[as.logical(1-deltas)]

# fix censoring time
# pct <- 0.75
# idx <- order(yy)[round(n*pct)]
# cc <- yy[idx]
# deltas <- yy<=cc
# sum(1-deltas)/n # percentage censored
# y <- yy
# y[as.logical(1-deltas)] <- cc

# gird plot
numpoints <- 100
mu.seq <- seq(-.5+mu0,.5+mu0,length.out = numpoints)
logel.seq <- rep(NA,numpoints)
for (ii in 1:numpoints) {
  G <- mr.evalG(y,X,mu.seq[ii])
  epsilons <- y - c(X * mu.seq[ii])
  logel.seq[ii] <- logEL(G,deltas,epsilons)
}
logelmode <- plotEL(mu.seq, logel.seq, mu0, NA, expression(mu))

nsamples <- 20000
nburn <- 3000
betaInit <- c(lm(y ~ 1)$coefficients) 
betaInit
sigs <- rep(0.3,1)
system.time(
  qrout <- mr_cens.post(y, X, deltas, nsamples, nburn, betaInit, sigs)
)

mu_chain <- qrout$beta_chain
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

# ---- 2-d problem (1 intercept, 1 slope) ----
n <- 50
p <- 2
X1 <- matrix(rnorm(n),n,1)
X <- cbind(1,X1)
eps <- rnorm(n) # N(0,1) error term
beta_I <- 1
beta_S <- 1.5
yy <- beta_I + c(X1 %*% beta_S) + eps 
beta0 <- c(beta_I, beta_S)
# plot(X1,y,cex=0.3)

# random censoring
cc <- rnorm(n,mean=1.5*abs(yy),sd=1)
deltas <- yy<=cc
y <- yy
sum(1-deltas)/n
y[as.logical(1-deltas)] <- cc[as.logical(1-deltas)]

# grid plot of conditionals: beta1|beta2 and beta2|beta1
numpoints <- 100
beta1.seq <- seq(beta0[1]-.5,beta0[1]+.5,length.out = numpoints)
beta2.seq <- seq(beta0[2]-.5,beta0[2]+.5,length.out = numpoints)
beta.seq <- cbind(beta1.seq,beta2.seq)
logel.seq <- matrix(rep(NA,2*numpoints),2,numpoints)
for (ii in 1:numpoints) {
  G <- mr.evalG(y,X,c(beta1.seq[ii],beta0[2]))
  epsilons <- y - c(X %*% c(beta1.seq[ii],beta0[2]))
  logel.seq[1,ii] <- logEL(G,deltas,epsilons)
  
  G <- mr.evalG(y,X,c(beta0[1],beta2.seq[ii]))
  epsilons <- y - c(X %*% c(beta0[1],beta2.seq[ii]))
  logel.seq[2,ii] <- logEL(G,deltas,epsilons)
}
logelmode1 <- plotEL(beta1.seq, logel.seq[1,], beta0[1], NA, expression(beta[0]))
logelmode2 <- plotEL(beta2.seq, logel.seq[2,], beta0[2], NA, expression(beta[1]))

# calculate marginal posterior
# Beta.seq <- as.matrix(expand.grid(beta1.seq, beta2.seq))
# logel.mat <- apply(Beta.seq, 1, function(bb) {
#   G <- mr.evalG(y,X,c(bb[1],bb[2]))
#   logEL(G)
# })
# logel.mat <- matrix(logel.mat, numpoints, numpoints)
# el.mat <- exp(logel.mat - max(logel.mat))
# logel.marg <- log(cbind(beta1 = rowSums(el.mat), beta2 = colSums(el.mat)))

nsamples <- 20000
nburn <- 3000
betaInit <- c(lm(y ~ X1)$coefficients)
betaInit
sigs <- rep(0.1,2)
# Note: For matching the conditional grid plot and histogram, pay attention to 
#   their inital values, the fixed params need to be the same.
RvDoMcmc <- c(1,0)
# k <- which(as.logical(RvDoMcmc))
# betaInit[-k] <- beta0[-k] # fix other params at their true values
system.time(
  qrout <- mr_cens.post(y, X, deltas, nsamples, nburn, c(betaInit[1],beta0[2]), sigs, RvDoMcmc)
)
beta_chain <- qrout$beta_chain
beta_paccept <- qrout$paccept
beta_paccept
plot(beta_chain[1,], xlab = expression(beta[0]), ylab = 'EL', type='l')
plot(beta_chain[2,], xlab = expression(beta[1]), ylab = 'EL', type='l')

# overlay marginal gird plot to histogram
# intercept
hist(beta_chain[1,],breaks=50,freq=FALSE,
     xlab = expression(beta[0]),main='')
# marginal line
lines(beta1.seq, norm_pdf(logel.marg[,1], beta1.seq),
      cex=0.1, col = 'red', type='l')
# conditional line
lines(beta1.seq, norm_pdf(logel.seq[1,], beta1.seq),
      cex=0.1, col = 'red', type='l')
abline(v=mean(beta_chain[1,]), col='blue')
legend('topright',legend=c(expression('grid plot'),
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
abline(v=mean(beta_chain[2,]), col='blue')
legend('topright',legend=c(expression('grid plot & mode'),
                           expression('sample mean')),
       lty = c(1,1), col = c('red','blue'), cex = 0.6)
