# check that logEL plot under location-scale models
require(bayesEL)
require(optimCheck)
source("hlm-functions.R")
source("../testthat/el-utils.R")

#---- mean reg: X and Z both dim 1----
# dimensions
n <- 100 # number of observations
numpoints <- 100

p <- 1
q <- 1
X <- matrix(rep(1,n),n,p)
# X <- matrix(rnorm(n),n,p)
# Z <- matrix(rep(1,n),n,q)
Z <- matrix(rnorm(n),n,q)

beta0 <- rnorm(p)
gamma0 <- rnorm(q)
# beta0 <- 2
# gamma0 <- -0.01
theta0 <- c(beta0,gamma0)
y <- c(X %*% beta0 + exp(Z %*% gamma0)*rnorm(n))

plot(y,cex=0.3)
beta.seq <- seq(-.5+beta0, .5+beta0, length.out = numpoints)
gamma.seq <- seq(-.5+gamma0, .5+gamma0, length.out = numpoints)
theta.seq <- rbind(beta.seq, gamma.seq)

logel.seq <- matrix(rep(NA,2*numpoints),2,numpoints)
for (ii in 1:numpoints) {
  G <- mrls.evalG(y,X,Z,beta.seq[ii],gamma0)
  logel.seq[1,ii] <- logEL(G)
  G <- mrls.evalG(y,X,Z,beta0,gamma.seq[ii])
  logel.seq[2,ii] <- logEL(G)
}
plot(beta.seq,logel.seq[1,], type='l')
abline(v=beta0,col='red')
plot(gamma.seq, logel.seq[2,], type='l')
abline(v=gamma0,col='red')

# solve by hlm: 
theta.hat <- hlm.fit(y = y, X = X, W = Z)
# quick check
rbind(true = beta0, est = theta.hat$beta) # beta
rbind(true = gamma0, est = theta.hat$gamma) # gamma

# optim check
loglik <- function(theta) {
  hlm.loglik(beta = theta[1:p], gamma = theta[p+1:q],
             y = y, X = X, W = Z)
}
theta.names <- c(paste0("beta[", 1:p-1, "]"), paste0("gamma[", 1:q-1, "]"))
theta.names <- parse(text = theta.names) # convert to greek symbols
optim_proj(fun = loglik,
           xsol = c(theta.hat$beta, theta.hat$gamma),
           xnames = theta.names)

# calculate marginal posterior
Theta.seq <- as.matrix(expand.grid(beta.seq, gamma.seq))
logel.mat <- apply(Theta.seq, 1, function(bb) {
  G <- mrls.evalG(y,X,Z,bb[1],bb[2])
  logEL(G)
})
logel.mat <- matrix(logel.mat, numpoints, numpoints)
el.mat <- exp(logel.mat - max(logel.mat))
logel.marg <- log(cbind(beta1 = rowSums(el.mat), beta2 = colSums(el.mat)))

# mcmc
nsamples <- 20000
nburn <- 3000
betaInit <- theta.hat$beta
gammaInit <- theta.hat$gamma
sigs <- c(0.25,0.23)

system.time(
  postout <- mrls.post(y,X,Z,nsamples,nburn,betaInit,gammaInit,sigs)
)
theta_chain <- postout$Theta_chain
theta_accept <- postout$paccept
theta_accept
hist(theta_chain[1,],breaks=50,freq=FALSE,
     xlab = expression(beta[0]), main='')
lines(beta.seq, norm_pdf(logel.marg[,1], beta.seq),
      cex=0.1, col = 'red', type='l')
hist(theta_chain[2,],breaks=50,freq=FALSE,
     xlab = expression(beta[1]), main='')
lines(gamma.seq, norm_pdf(logel.marg[,2], beta.seq),
      cex=0.1, col = 'red', type='l')

# ---- mean reg: X and Z both dim 2 ---- 
# dimensions
n <- 100 # number of observations
numpoints <- 100

p <- 2
q <- 2
X <- cbind(1,rnorm(n))
# X <- matrix(rnorm(p*n),n,p)
Z <- matrix(rnorm(q*n),n,q)
# Z <- cbind(1,rnorm(n))

# beta0 <- rnorm(p)
# gamma0 <- rnorm(q)
beta0 <- c(1,2)
gamma0 <- c(-0.8, -0.2)
y <- c(X %*% beta0 + exp(Z %*% gamma0)*rnorm(n))
plot(y,cex=0.3)
beta1.seq <- seq(-.5+beta0[1], .5+beta0[1], length.out = numpoints)
beta2.seq <- seq(-.5+beta0[2], .5+beta0[2], length.out = numpoints)
gamma1.seq <- seq(-.5+gamma0[1], .5+gamma0[1], length.out = numpoints)
gamma2.seq <- seq(-.5+gamma0[2], .5+gamma0[2], length.out = numpoints)

logel.seq <- matrix(rep(NA,p*q*numpoints),p*q,numpoints)
for (ii in 1:numpoints) {
  G <- mrls.evalG(y,X,Z,c(beta1.seq[ii],beta0[2]),gamma0)
  logel.seq[1,ii] <- logEL(G)
  G <- mrls.evalG(y,X,Z,c(beta0[1],beta2.seq[ii]),gamma0)
  logel.seq[2,ii] <- logEL(G)
  G <- mrls.evalG(y,X,Z,beta0,c(gamma1.seq[ii],gamma0[2]))
  logel.seq[3,ii] <- logEL(G)
  G <- mrls.evalG(y,X,Z,beta0,c(gamma0[1],gamma2.seq[ii]))
  logel.seq[4,ii] <- logEL(G)
}
plot(beta1.seq,logel.seq[1,], type='l')
abline(v=beta0[1],col='red')
plot(beta2.seq,logel.seq[2,], type='l')
abline(v=beta0[2],col='red')
plot(gamma1.seq, logel.seq[3,], type='l')
abline(v=gamma0[1],col='red')
plot(gamma2.seq, logel.seq[4,], type='l')
abline(v=gamma0[2],col='red')

# solve by hlm
theta.hat <- hlm.fit(y = y, X = X, W = Z)
# quick check
rbind(true = beta0, est = theta.hat$beta) # beta
rbind(true = gamma0, est = theta.hat$gamma) # gamma

loglik <- function(theta) {
  hlm.loglik(beta = theta[1:p], gamma = theta[p+1:q],
             y = y, X = X, W = Z)
}

# optim check
theta.names <- c(paste0("beta[", 1:p-1, "]"), paste0("gamma[", 1:q-1, "]"))
theta.names <- parse(text = theta.names) # convert to greek symbols
optim_proj(fun = loglik,
           xsol = c(theta.hat$beta, theta.hat$gamma),
           xnames = theta.names)

# mcmc
nsamples <- 20000
nburn <- 3000
betaInit <- theta.hat$beta
gammaInit <- theta.hat$gamma
sigs <- c(0.15,0.1,0.23,0.23)

system.time(
  postout <- mrls.post(y,X,Z,nsamples,nburn,betaInit,gammaInit,sigs)
)
theta_chain <- postout$Theta_chain
theta_accept <- postout$paccept
theta_accept
hist(theta_chain[1,],breaks=50,freq=FALSE,
     xlab = expression(beta[0]), main='')
hist(theta_chain[2,],breaks=50,freq=FALSE,
     xlab = expression(beta[1]), main='')
hist(theta_chain[3,],breaks=50,freq=FALSE,
     xlab = expression(gamma[0]), main='')
hist(theta_chain[4,],breaks=50,freq=FALSE,
     xlab = expression(gamma[1]), main='')

# --------------------------------------------------------------------------- 
# test qr.evalG and qrls.evalG for multiple quantile levels
n <- 20
alphas <- c(0.25,0.5)
Beta <- matrix(c(0.5,1,1,1),2,2)
x <- rnorm(n)
y <- 1 + x + rnorm(n)
X <- cbind(1,x)
qr.evalG(y,X,alphas,Beta)
qr.evalG(y,X,alphas[1],Beta[,1])
qr.evalG(y,X,alphas[2],Beta[,2])

# Beta <- cbind(rnorm(p),rnorm(p))
# Gamma <- cbind(rnorm(q),rnorm(q))
# t(qrls.evalG(y,X,Z,alphas,Beta,Gamma))
# --------------------------------------------------------------------------- 

# ---- quant reg: X and Z both dim 1 with 1 quant level ----
# dimensions
n <- 100 # number of observations
numpoints <- 100

p <- 1
q <- 1
# X <- matrix(rep(1,n),n,p)
X <- matrix(rnorm(n),n,p)
# Z <- matrix(rep(1,n),n,q)
Z <- matrix(rnorm(n),n,q)
# Z <- X

# beta0 <- rnorm(p)
# beta0
# gamma0 <- rnorm(q)
# gamma0
# beta0 <- 0.2910
# gamma0 <- -0.5644
beta0 <- 1
gamma0 <- -0.5
alpha <- 0.5
eps <- rnorm(n)+0.1
quantile(eps,alpha)
nu0 <- 0.1
y <- c(X %*% beta0 + exp(Z %*% gamma0)*eps)
plot(y~X, cex=0.3)

# solve by hlm: 
theta.hat <- hlm.fit(y = y, X = X, W = Z)
# quick check
rbind(true = beta0, est = theta.hat$beta) # beta
rbind(true = gamma0, est = theta.hat$gamma) # gamma

eps_new <- c((y - X %*% theta.hat$beta)/exp(Z %*% theta.hat$gamma))

beta.seq <- seq(-.5+beta0, .5+beta0, length.out = numpoints)
gamma.seq <- seq(-.5+gamma0, .5+gamma0, length.out = numpoints)
nu.seq <- seq(-.5+nu0, .5+nu0, length.out = numpoints)
theta.seq <- rbind(beta.seq, gamma.seq)

logel.seq <- matrix(rep(NA,3*numpoints),3,numpoints)
for (ii in 1:numpoints) {
  G <- qrls.evalG(y,X,Z,alpha,beta.seq[ii],gamma0,nu0)
  logel.seq[1,ii] <- logEL(G)
  G <- qrls.evalG(y,X,Z,alpha,beta0,gamma.seq[ii],nu0)
  logel.seq[2,ii] <- logEL(G)
  G <- qrls.evalG(y,X,Z,alpha,beta0,gamma0,nu.seq[ii])
  logel.seq[3,ii] <- logEL(G)
}
plot(beta.seq,logel.seq[1,], type='l')
abline(v=beta0,col='red')
plot(gamma.seq, logel.seq[2,], type='l')
abline(v=gamma0,col='red')
plot(nu.seq, logel.seq[3,], type='l')
abline(v=nu0,col='red')

# 3-way
Theta.seq <- as.matrix(expand.grid(beta.seq, gamma.seq, nu.seq))
logel.mat <- apply(Theta.seq, 1, function(bb) {
  G <- qrls.evalG(y,X,Z,alpha,bb[1],bb[2],bb[3])
  logEL(G)
})
logel.mat <- array(logel.mat, c(numpoints, numpoints,numpoints))
el.mat <- exp(logel.mat - max(logel.mat))
logel.marg1 <- log(apply(el.mat, MARGIN=1, sum))
logel.marg2 <- log(apply(el.mat, MARGIN=2, sum))
logel.marg3 <- log(apply(el.mat, MARGIN=3, sum))

# mcmc
nsamples <- 20000
nburn <- 3000
betaInit <- theta.hat$beta
gammaInit <- theta.hat$gamma
nuInit <- quantile(eps_new,alpha)
# sigs <- c(0.45,0.3,1)
sigs <- c(0.25,0.32,0.26)

system.time(
  postout <- qrls.post(y,X,Z,alpha,nsamples,nburn,betaInit,gammaInit,nuInit,sigs)
)
theta_chain <- postout$Theta_chain
theta_accept <- postout$paccept
theta_accept
hist(theta_chain[1,],breaks=50,freq=FALSE,
     xlab = expression(beta), main='')
lines(beta.seq, norm_pdf(logel.marg1, beta.seq),
      cex=0.1, col = 'red', type='l')
abline(v=beta0, col='red')
abline(v=mean(theta_chain[1,]), col='blue')
legend('topright',legend=c(expression('grid plot & true param'),
                           expression('sample mean')),
       lty = c(1,1), col = c('red','blue'), cex = 0.6)
hist(theta_chain[2,],breaks=50,freq=FALSE,
     xlab = expression(gamma), main='')
lines(gamma.seq, norm_pdf(logel.marg2, gamma.seq),
      cex=0.1, col = 'red', type='l')
abline(v=gamma0, col='red')
abline(v=mean(theta_chain[2,]), col='blue')
legend('topright',legend=c(expression('grid plot & true param'),
                           expression('sample mean')),
       lty = c(1,1), col = c('red','blue'), cex = 0.6)
hist(theta_chain[3,],breaks=50,freq=FALSE,
     xlab = expression(nu), main='')
lines(nu.seq, norm_pdf(logel.marg3, nu.seq),
      cex=0.1, col = 'red', type='l')
abline(v=nu0, col='red')
abline(v=mean(theta_chain[3,]), col='blue')
legend('topright',legend=c(expression('grid plot & true param'),
                           expression('sample mean')),
       lty = c(1,1), col = c('red','blue'), cex = 0.6)

# fit with only location model
XX <- cbind(1,X)
betaInit <- c(rq(y~X, tau=alpha)$coefficients)
betaInit
system.time(
  postout_l <- qr.post(y,XX,alpha,nsamples,nburn,betaInit,c(0.2,0.2))
)
hist(postout_l$Beta_chain[1,],breaks=50,freq=FALSE,
     xlab = expression(beta[0]), main='')
abline(v=nu0, col='red')
abline(v=mean(postout_l$Beta_chain[1,]), col='blue')
legend('topright',legend=c(expression('true param'),
                           expression('sample mean')),
       lty = c(1,1), col = c('red','blue'), cex = 0.6)
hist(postout_l$Beta_chain[2,],breaks=50,freq=FALSE,
     xlab = expression(beta[1]), main='')
abline(v=nu0, col='red')
abline(v=mean(postout_l$Beta_chain[2,]), col='blue')
legend('topright',legend=c(expression('true param'),
                           expression('sample mean')),
       lty = c(1,1), col = c('red','blue'), cex = 0.6)

# repeat the above and calculate the MSE of the slope parameter
# alpha <- c(0.25,0.5,0.75)
# alpha <- c(0.9,0.925,0.95)
alpha <- 0.75
nrep <- 120
n <- 200
nsamples <- 15000
nburn <- 10000
sigs <- c(0.2,0.3,0.25)

beta0 <- 1
gamma0 <- -0.5
nu0 <- 0.1+qnorm(alpha)

beta_true <- c(nu0,beta0)

RQest <- matrix(rep(NA,2*nrep),2,nrep)
ELest <- matrix(rep(NA,2*nrep),2,nrep)
system.time(
  for (ii in 1:nrep) {
    eps <- rnorm(n)+0.1
    X <- matrix(rnorm(n),n,1)
    Z <- X
    y <- c(X %*% beta0 + exp(Z %*% gamma0)*eps)
    theta.hat <- tryCatch(hlm.fit(y = y, X = X, W = Z), error = function() next)
    # theta.hat <- hlm.fit(y = y, X = X, W = Z)
    betaInit <- theta.hat$beta
    gammaInit <- theta.hat$gamma
    eps_new <- c((y - X %*% theta.hat$beta)/exp(Z %*% theta.hat$gamma))
    nuInit <- quantile(eps_new,alpha)
    postout <- qrls.post(y,X,Z,alpha,nsamples,nburn,betaInit,gammaInit,nuInit,sigs)
    ELest[,ii] <- rowMeans(postout$Theta_chain[c(3,1),])
    XX <- cbind(1,X)
    betaInit <- c(rq(y~X, tau=alpha)$coefficients)
    postout_l <- qr.post_adapt(y,XX,alpha,nsamples,nburn,betaInit,c(0.2,0.2))
    RQest[,ii] <- rowMeans(postout_l$beta_chain[1:2,])
    if (ii %% 5 == 0) {
      message("ii = ", ii)
      message("cumu MSE of ELest = ", sum((ELest[,1:ii] - beta_true)^2)/ii)
      message("cumu MSE of RQest = ", sum((RQest[,1:ii] - beta_true)^2)/ii)
    }
  }
)

sum((ELest[1,1:100] - beta_true)^2)/nrep
sum((ELest[2,1:100] - beta_true)^2)/nrep
sum((RQest[1,1:100] - beta_true)^2)/nrep
sum((RQest[2,1:100] - beta_true)^2)/nrep


# TODO: not done yet
# ---- quant reg: X and Z both dim 1 with 2 quant level ----
library(quantreg)
# dimensions
n <- 200 # number of observations
numpoints <- 100

p <- 1
q <- 1
# X <- matrix(rep(1,n),n,p)
X <- matrix(rnorm(n),n,p)
# Z <- matrix(rep(1,n),n,q)
Z <- matrix(rnorm(n),n,q)

# beta0 <- rnorm(p)
# gamma0 <- rnorm(q)
beta0 <- 2
gamma0 <- -0.01
theta0 <- c(beta0,gamma0)
eps <- rnorm(n)
y <- c(X %*% beta0 + exp(Z %*% gamma0)*eps)
plot(y~X, cex=0.3)
alpha <- c(0.25,0.75)

theta.hat <- hlm.fit(y,X,Z)
betaInit <- theta.hat$beta
gammaInit <- theta.hat$gamma
BetaInit <- cbind(betaInit,betaInit)
GammaInit <- cbind(gammaInit,gammaInit)
eps_new <- c((y - X %*% theta.hat$beta)/exp(Z %*% theta.hat$gamma))
nuInit <- c(quantile(eps_new,alpha[1]),quantile(eps_new,alpha[2]))
NuInit <- nuInit
Sigs <- rbind(c(0.1,0.1),c(0.1,0.15),c(0.1,0.15))
nsamples <- 20000
nburn <- 3000

G <- qrls.evalG(y,X,Z,alpha,BetaInit,GammaInit,NuInit)
G1 <- qrls.evalG(y,X,Z,alpha[1],BetaInit[,1],GammaInit[,1],NuInit[1])
G2 <- qrls.evalG(y,X,Z,alpha[2],BetaInit[,2],GammaInit[,2],NuInit[2])
lambdaNR(G)
lambdaNR(G1)
lambdaNR(G2)
postout1 <- qrls.post(y,X,Z,alpha[1],nsamples,nburn,
                     BetaInit[,1],GammaInit[,1],NuInit[1],Sigs[,1])
postout2 <- qrls.post(y,X,Z,alpha[2],nsamples,nburn,
                      BetaInit[,2],GammaInit[,2],NuInit[2],Sigs[,2])
hist(postout1$Theta_chain[1,],breaks=50,freq=FALSE,
     xlab = expression(beta[0]), main='')
hist(postout1$Theta_chain[2,],breaks=50,freq=FALSE,
     xlab = expression(gamma[0]), main='')
hist(postout2$Theta_chain[1,],breaks=50,freq=FALSE,
     xlab = expression(beta[1]), main='')
hist(postout2$Theta_chain[2,],breaks=50,freq=FALSE,
     xlab = expression(gamma[1]), main='')

system.time(
  postout <- qrls.post(y,X,Z,alpha,nsamples,nburn,
                       BetaInit,GammaInit,NuInit,Sigs)
)
# NuInit = -0.7763309  0.6910551 
theta_chain <- postout$Theta_chain
paccept <- postout$paccept
paccept
hist(theta_chain[1,],breaks=50,freq=FALSE,
     xlab = expression(beta[0]), main='')
hist(theta_chain[2,],breaks=50,freq=FALSE,
     xlab = expression(gamma[0]), main='')
hist(theta_chain[4,],breaks=50,freq=FALSE,
     xlab = expression(beta[1]), main='')
hist(theta_chain[5,],breaks=50,freq=FALSE,
     xlab = expression(gamma[1]), main='')

# --------------------------------------------------------------------------- 
# ---- 1-d problem with 2 quantile levels ---- 
n <- 300
mu <- 1
alphas <- c(0.9,0.925,0.95)
X <- matrix(rep(1,n), n, 1) # each row of X is one observation
eps <- rnorm(n) # N(0,1) error term
quantile(eps, alphas[1])
quantile(eps, alphas[2])
y <- c(X * mu) + eps
mu01 <- mu + qnorm(alphas[1], lower.tail = TRUE) # true param assuming alpha-tail of eps is 0
mu01
mu02 <- mu + qnorm(alphas[2], lower.tail = TRUE) # true param assuming alpha-tail of eps is 0
mu02
mu03 <- mu + qnorm(alphas[3], lower.tail = TRUE) # true param assuming alpha-tail of eps is 0
mu03

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
  # logel.seq2[ii] <- logEL(G = G)
  logel.seq2[ii] <- logEL_shrink(G = G, mu01, mu.seq2[ii], 0.5)
}
logelmode2 <- plotEL(mu.seq2, logel.seq2, mu02, quantile(y,alphas[2]), expression(mu))
logelmode2

numpoints <- 100
mu.seq3 <- seq(-.5+mu03,.5+mu03,length.out = numpoints)
logel.seq3 <- rep(NA,numpoints)
for (ii in 1:numpoints) {
  G <- qr.evalG(y,X,alphas[3],mu.seq3[ii])
  # logel.seq2[ii] <- logEL(G = G)
  logel.seq3[ii] <- logEL_shrink(G = G, mu01, mu.seq2[ii], 0.5)
}
logelmode2 <- plotEL(mu.seq3, logel.seq3, mu03, quantile(y,alphas[3]), expression(mu))
logelmode2

nsamples <- 20000
nburn <- 3000
library(quantreg)
betaInit1 <- c(rq(y ~ 1, tau = alphas[1], method = 'fn')$coefficients)
betaInit1
betaInit2 <- c(rq(y ~ 1, tau = alphas[2], method = 'fn')$coefficients)
betaInit2
betaInit3 <- c(rq(y ~ 1, tau = alphas[3], method = 'fn')$coefficients)
betaInit3
BetaInit <- cbind(betaInit1,betaInit2,betaInit3)
# BetaInit <- cbind(mu01,mu02)
sigs1 <- rep(0.15,1)
sigs2 <- rep(0.15,1)
sigs3 <- rep(0.15,1)
Sigs <- cbind(sigs1,sigs2,sigs3)
qrout <- qr.post_shrink(y, X, alphas[2:3], nsamples, nburn, rbind(BetaInit[,2:3]), rbind(Sigs[,2:3]), 1)
# qrout <- qr.post(y, X, alphas, nsamples, nburn, BetaInit, Sigs)

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
abline(v=mu01, col='red')
abline(v=mean(mu_chain[1,]), col='blue')
legend('topright',legend=c(expression('grid plot & true param'),
                           expression('sample mean')),
       lty = c(1,1), col = c('red','blue'), cex = 0.6)

hist(mu_chain[2,],breaks=50,freq=FALSE,
     xlab = expression(mu), main='')
lines(mu.seq2, norm_pdf(logel.seq2, mu.seq2),
      cex=0.1, col = 'red', type='l')
abline(v=mu02, col='red')
abline(v=mean(mu_chain[2,]), col='blue')
legend('topright',legend=c(expression('grid plot & true param'),
                           expression('sample mean')),
       lty = c(1,1), col = c('red','blue'), cex = 0.6)

hist(mu_chain[3,],breaks=50,freq=FALSE,
     xlab = expression(mu), main='')
lines(mu.seq3, norm_pdf(logel.seq3, mu.seq3),
      cex=0.1, col = 'red', type='l')
abline(v=mu03, col='red')
abline(v=mean(mu_chain[3,]), col='blue')
legend('topright',legend=c(expression('grid plot & true param'),
                           expression('sample mean')),
       lty = c(1,1), col = c('red','blue'), cex = 0.6)

