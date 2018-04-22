require(bayesEL)
require(optimCheck)
source("hlm-functions.R")
source("../testthat/el-utils.R")
# ---- quant reg: X and Z both dim 1 with 1 quant level ----
# dimensions
n <- 100 # number of observations
p <- 1
q <- 1
# X <- matrix(rep(1,n),n,p)
X <- matrix(rnorm(n),n,p)
# Z <- matrix(rep(1,n),n,q) # Z should not contain intercept too ???
Z <- matrix(rnorm(n),n,q)
# Z <- X

beta0 <- rnorm(p)
gamma0 <- rnorm(q)
# beta0 <- 1
# gamma0 <- -0.5
alpha <- 0.75
eps <- rnorm(n)
quantile(eps,alpha)
nu0 <- qnorm(alpha)
y <- c(X %*% beta0 + exp(Z %*% gamma0)*eps)
plot(y~X, cex=0.3)

numpoints <- 50
beta.seq <- seq(-.5+beta0, .5+beta0, length.out = numpoints)
gamma.seq <- seq(-.5+gamma0, .5+gamma0, length.out = numpoints)
nu.seq <- seq(-.5+nu0, .5+nu0, length.out = numpoints)
theta.seq <- rbind(beta.seq, gamma.seq)

# calculate the marginal posterior distribution
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
theta.hat <- hlm.fit(y = y, X = X, W = Z) # solve by hlm
betaInit <- theta.hat$beta
gammaInit <- theta.hat$gamma
eps_new <- c((y - X %*% theta.hat$beta)/exp(Z %*% theta.hat$gamma))
nuInit <- quantile(eps_new,alpha)
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
