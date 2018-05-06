require(bayesEL)
source("hlm-functions.R")
source("../testthat/el-utils.R")

##### 1 quantile level #####

# ---- X and Z both dim 1 (works fine for all eps here) ----
# dimensions
n <- 200 # number of observations
p <- 1
q <- 1

# covariates
# X <- matrix(rep(1,n),n,p)
X <- matrix(rnorm(n),n,p)
# Z <- matrix(rep(1,n),n,q) # Z should not contain intercept
Z <- matrix(rnorm(n),n,q)
# Z <- X

# parameters
beta0 <- rnorm(p)
beta0
gamma0 <- rnorm(q)
gamma0
# beta0 <- 0.5
# gamma0 <- -0.5
alpha <- 0.75

# normal(0,1) error
eps <- rnorm(n)
nu0 <- qnorm(alpha)

# t error with mean 0 var 1
df <- 5
v <- df/(df-2)
eps <- rt(n, df=df)/sqrt(v)
nu0 <- qt(alpha,df=df)/sqrt(v)

# chi-sqr with mean 0 var 1
df <- 3
m <- df
v <- 2*df
eps <- (rchisq(n, df=df)-m)/sqrt(v)
nu0 <- (qchisq(alpha, df)-m)/sqrt(v)

# log-normal with mean 0 var 1
# {\displaystyle [\exp(\sigma ^{2})-1]\exp(2\mu +\sigma ^{2})}
mn <- 0
sn <- 1
m <- exp(mn+sn^2/2)
v <- (exp(sn^2)-1)*exp(2*mn+sn^2)
eps <- (rlnorm(n,mn,sn)-m)/sqrt(v)
nu0 <- (qlnorm(alpha,mn,sn)-m)/sqrt(v)

# response
y <- c(X %*% beta0 + exp(0.5 * Z %*% gamma0)*eps)
plot(y~X, cex=0.3)

# calculate the marginal posterior distributions
numpoints <- 50
beta.seq <- seq(-1+beta0, 1+beta0, length.out = numpoints)
gamma.seq <- seq(-1+gamma0, 1+gamma0, length.out = numpoints)
nu.seq <- seq(-1+nu0, 1+nu0, length.out = numpoints)

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

# calculate the conditional posterior distributions
numpoints <- 100
beta.seq <- seq(-1+beta0,1+beta0,length.out = numpoints)
logel.seq1 <- rep(NA,numpoints)
for (ii in 1:numpoints) {
  G <- qrls.evalG(y,X,Z,alpha,beta.seq[ii],gamma0/2,nu0)
  logel.seq1[ii] <- logEL(G = G)
}
logelmode1 <- plotEL(beta.seq, logel.seq1, beta0, NA, expression(beta[0]))

gamma.seq <- seq(-1+gamma0/2,1+gamma0/2,length.out = numpoints)
logel.seq2 <- rep(NA,numpoints)
for (ii in 1:numpoints) {
  G <- qrls.evalG(y,X,Z,alpha,beta0,gamma.seq[ii],nu0)
  logel.seq2[ii] <- logEL(G = G)
}
logelmode2 <- plotEL(gamma.seq, logel.seq2, gamma0/2, NA, expression(gamma[0]))

nu.seq <- seq(-1+nu0,1+nu0,length.out = numpoints)
logel.seq3 <- rep(NA,numpoints)
for (ii in 1:numpoints) {
  G <- qrls.evalG(y,X,Z,alpha,beta0,gamma0/2,nu.seq[ii])
  logel.seq3[ii] <- logEL(G = G)
}
logelmode3 <- plotEL(nu.seq, logel.seq3, nu0, NA, expression(nu))

# mcmc
nsamples <- 20000
nburn <- 3000
theta.hat <- hlm.fit(y = y, X = X, W = Z) # solve by hlm
betaInit <- theta.hat$beta
gammaInit <- theta.hat$gamma/2
eps_new <- c((y - X %*% theta.hat$beta)*exp(-0.5*Z %*% theta.hat$gamma))
nuInit <- quantile(eps_new,alpha)
sigs <- c(0.2,0.2,0.2)
RvDoMcmc <- c(0,1,0)
system.time(
  # postout <- qrls.post(y,X,Z,alpha,nsamples,nburn,beta0,gamma0/2,nu0,sigs,RvDoMcmc)
  postout <- qrls.post_adapt(y,X,Z,alpha,nsamples,nburn,beta0,gammaInit,nu0,sigs,RvDoMcmc)
)
theta_chain <- postout$theta_chain
theta_accept <- postout$paccept
theta_accept

# mixing of the chain
plot(theta_chain[2,],type = 'l')

# histograms
hist(theta_chain[1,],breaks=50,freq=FALSE,
     xlab = expression(beta), main='')
# marginal line
lines(beta.seq, norm_pdf(logel.marg1, beta.seq),
      cex=0.1, col = 'red', type='l')
# conditional line
lines(beta.seq, norm_pdf(logel.seq1, beta.seq),
      cex=0.1, col = 'blue', type='l')
abline(v=beta0, col='red')
abline(v=mean(theta_chain[1,]), col='blue')
legend('topright',legend=c(expression('true param'),
                           expression('sample mean')),
       lty = c(1,1), col = c('red','blue'), cex = 0.6)

hist(theta_chain[2,],breaks=50,freq=FALSE,
     xlab = expression(gamma), main='')
# marginal line
lines(gamma.seq, norm_pdf(logel.marg2, gamma.seq),
      cex=0.1, col = 'red', type='l')
# conditional line
lines(gamma.seq, norm_pdf(logel.seq2, gamma.seq),
      cex=0.1, col = 'blue', type='l')
abline(v=gamma0/2, col='red')
abline(v=mean(theta_chain[2,]), col='blue')
legend('topright',legend=c(expression('true param'),
                           expression('sample mean')),
       lty = c(1,1), col = c('red','blue'), cex = 0.6)

hist(theta_chain[3,],breaks=50,freq=FALSE,
     xlab = expression(nu), main='')
# marginal line
lines(nu.seq, norm_pdf(logel.marg3, nu.seq),
      cex=0.1, col = 'red', type='l')
# conditional line
lines(nu.seq, norm_pdf(logel.seq3, nu.seq),
      cex=0.1, col = 'blue', type='l')
abline(v=nu0, col='red')
abline(v=mean(theta_chain[3,]), col='blue')
legend('topright',legend=c(expression('true param'),
                           expression('sample mean')),
       lty = c(1,1), col = c('red','blue'), cex = 0.6)

# ---- X and Z both dim 2 ----
# dimensions
n <- 200 # number of observations
p <- 2
q <- 2
X <- cbind(rep(1,n),rnorm(n))
# X <- matrix(rnorm(p*n),n,p)
# Z <- matrix(rep(1,n),n,q) # Z should not contain an intercept
Z <- matrix(rnorm(q*n),n,q)

# parameters
beta0 <- rnorm(p)
beta0
gamma0 <- rnorm(q)
gamma0
alpha <- 0.75

# normal(0,1) error
eps <- rnorm(n)
nu0 <- qnorm(alpha)

# t error with mean 0 var 1
df <- 5
v <- df/(df-2)
eps <- rt(n, df=df)/sqrt(v)
nu0 <- qt(alpha,df=df)/sqrt(v)

# chi-sqr with mean 0 var 1
df <- 3
m <- df
v <- 2*df
eps <- (rchisq(n, df=df)-m)/sqrt(v)
nu0 <- (qchisq(alpha, df)-m)/sqrt(v)

# log-normal with mean 0 var 1
# {\displaystyle [\exp(\sigma ^{2})-1]\exp(2\mu +\sigma ^{2})}
mn <- 0
sn <- 1
m <- exp(mn+sn^2/2)
v <- (exp(sn^2)-1)*exp(2*mn+sn^2)
eps <- (rlnorm(n,mn,sn)-m)/sqrt(v)
nu0 <- (qlnorm(alpha,mn,sn)-m)/sqrt(v)

# response
y <- c(X %*% beta0 + exp(0.5 * Z %*% gamma0)*eps)
# plot(y~X[,2], cex=0.3)

# for plotting conditional curves
numpoints <- 100

beta.seq1 <- seq(-1+beta0[1],1+beta0[1],length.out = numpoints)
logel.seq1 <- rep(NA,numpoints)
for (ii in 1:numpoints) {
  G <- qrls.evalG(y,X,Z,alpha,c(beta.seq1[ii],beta0[2]),gamma0/2,nu0)
  logel.seq1[ii] <- logEL(G = G)
}
logelmode1 <- plotEL(beta.seq1, logel.seq1, beta0[1], NA, expression(beta[0]))

beta.seq2 <- seq(-1+beta0[2],1+beta0[2],length.out = numpoints)
logel.seq2 <- rep(NA,numpoints)
for (ii in 1:numpoints) {
  G <- qrls.evalG(y,X,Z,alpha,c(beta0[1],beta.seq2[ii]),gamma0/2,nu0)
  logel.seq2[ii] <- logEL(G = G)
}
logelmode2 <- plotEL(beta.seq2, logel.seq2, beta0[2], NA, expression(beta[1]))

gamma.seq1 <- seq(-1+gamma0[1]/2,1+gamma0[1]/2,length.out = numpoints)
logel.seq3 <- rep(NA,numpoints)
for (ii in 1:numpoints) {
  G <- qrls.evalG(y,X,Z,alpha,beta0,c(gamma.seq1[ii],gamma0[2]/2),nu0)
  logel.seq3[ii] <- logEL(G = G)
}
logelmode3 <- plotEL(gamma.seq1, logel.seq3, gamma0[1]/2, NA, expression(gamma[0]))

gamma.seq2 <- seq(-1+gamma0[2]/2,1+gamma0[2]/2,length.out = numpoints)
logel.seq4 <- rep(NA,numpoints)
for (ii in 1:numpoints) {
  G <- qrls.evalG(y,X,Z,alpha,beta0,c(gamma0[1]/2,gamma.seq2[ii]),nu0)
  logel.seq4[ii] <- logEL(G = G)
}
logelmode4 <- plotEL(gamma.seq2, logel.seq4, gamma0[2]/2, NA, expression(gamma[1]))

nu.seq <- seq(-1+nu0,1+nu0,length.out = numpoints)
logel.seq5 <- rep(NA,numpoints)
for (ii in 1:numpoints) {
  G <- qrls.evalG(y,X,Z,alpha,beta0,gamma0/2,nu.seq[ii])
  logel.seq5[ii] <- logEL(G = G)
}
logelmode5 <- plotEL(nu.seq, logel.seq5, nu0, NA, expression(nu))

# mcmc
nsamples <- 20000
nburn <- 3000
# solve by hlm
theta.hat <- hlm.fit(y = y, X = X, W = cbind(1,Z))
betaInit <- theta.hat$beta
gammaInit <- theta.hat$gamma[2:(q+1)]/2
nuInit <- quantile((y-X %*% betaInit)*exp(-0.5*Z %*% gammaInit),alpha)
mwgSd <- c(0.2,0.2,0.2,0.2,0.2)

# choose which parameters to update, and others are fixed
RvDoMcmc <- rep(0,5)
RvDoMcmc[5] <- 1

system.time(
  # postout <- qrls.post(y,X,Z,alpha,nsamples,nburn,beta0,gamma0/2,nu0,mwgSd,RvDoMcmc)
  postout <- qrls.post_adapt(y,X,Z,alpha,nsamples,nburn,beta0,gamma0/2,nuInit,mwgSd,RvDoMcmc)
)
theta_chain <- postout$theta_chain
theta_accept <- postout$paccept
theta_accept

# mixing of the chains
plot(theta_chain[5,],type = 'l')

# histograms and conditional loglikelihood overlays
hist(theta_chain[1,],breaks=50,freq=FALSE,
     xlab = expression(beta[1]), main='')
# conditional line
lines(beta.seq1, norm_pdf(logel.seq1, beta.seq1),
      cex=0.1, col = 'blue', type='l')
abline(v=beta0[1],col='red')
abline(v=mean(theta_chain[1,]),col='blue')

hist(theta_chain[2,],breaks=50,freq=FALSE,
     xlab = expression(beta[2]), main='')
# conditional line
lines(beta.seq2, norm_pdf(logel.seq2, beta.seq2),
      cex=0.1, col = 'blue', type='l')
abline(v=beta0[2],col='red')
abline(v=mean(theta_chain[2,]),col='blue')

hist(theta_chain[3,],breaks=50,freq=FALSE,
     xlab = expression(gamma[1]), main='')
# conditional line
lines(gamma.seq1, norm_pdf(logel.seq3, gamma.seq1),
      cex=0.1, col = 'blue', type='l')
abline(v=gamma0[1]/2,col='red')
abline(v=mean(theta_chain[3,]),col='blue')

hist(theta_chain[4,],breaks=50,freq=FALSE,
     xlab = expression(gamma[2]), main='')
# conditional line
lines(gamma.seq2, norm_pdf(logel.seq4, gamma.seq2),
      cex=0.1, col = 'blue', type='l')
abline(v=gamma0[2],col='red')
abline(v=mean(theta_chain[4,]),col='blue')

hist(theta_chain[5,],breaks=50,freq=FALSE,
     xlab = expression(nu), main='')
# conditional line
lines(nu.seq, norm_pdf(logel.seq5, nu.seq),
      cex=0.1, col = 'blue', type='l')
abline(v=nu0,col='red')
abline(v=mean(theta_chain[5,]),col='blue')


