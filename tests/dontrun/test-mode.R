# import the helper functions
source("../testthat/el-utils.R")
source("../testthat/el-rfuns.R")
source("../testthat/el-model.R")
source("gen_eps.R")
source("mode-functions.R")

## mr ##

# 1-d problem 
n <- 100
mu0 <- 1
X <- matrix(rep(1,n), n, 1) # each row of X is one observation
eps <- rnorm(n) # N(0,1) error term
y <- c(X * mu0) + eps

mr.neglogEL_R(y,X,mu0)

numpoints <- 100
mu.seq <- seq(-.5+mu0,.5+mu0,length.out = numpoints)
logel.seq <- rep(NA,numpoints)
grad.seq <- rep(NA,numpoints)
for (ii in 1:numpoints) {
  temp <- -mr.neglogEL_R(y,X,mu.seq[ii])
  logel.seq[ii] <- temp
  grad.seq[ii] <- attributes(temp)$gradient
}
logelmode <- plotEL(mu.seq, logel.seq, mu0, mean(y), expression(mu))
plot(grad.seq,type='l')
abline(h=0,col='blue')

nlm(mr.neglogEL_R,1.2,y=y,X=X)

# 2-d problem
n <- 500
p <- 2
X0 <- matrix(rep(1,n),n,1)
X1 <- matrix(rnorm(n),n,1)
X <- cbind(X0,X1)
# eps <- rnorm(n) # N(0,1) error term
eps <- gen_eps(n,dist="nct",df=3,ncp=-1)
beta_intercept <- 1
beta_slope <- 1.5
y <- 1 + c(X1 %*% beta_slope) + eps 
beta0 <- c(beta_intercept, beta_slope)
beta.hat <- coef(lm(y~X-1))
beta.hat
mr.neglogEL_R(y,X,beta.hat)
nlmout <- nlm(mr.neglogEL_R,beta0,y=y,X=X,hessian=TRUE)
nlmout

# grid plot (conditionals: beta1|beta2 and beta2|beta1)
numpoints <- 100
beta1.seq <- seq(beta0[1]-1,beta0[1]+1,length.out = numpoints)
beta2.seq <- seq(beta0[2]-1,beta0[2]+1,length.out = numpoints)
# beta2.seq <- seq(1.45,1.55,length.out = numpoints)
logel.seq <- matrix(rep(NA,2*numpoints),2,numpoints)
for (ii in 1:numpoints) {
  if (ii %% 20 == 0) message("ii = ", ii)
  logel.seq[1,ii] <- -mr.neglogEL_R(y,X,c(beta1.seq[ii],beta0[2]))
  logel.seq[2,ii] <- -mr.neglogEL_R(y,X,c(beta0[1],beta2.seq[ii]))
}
logelmode1 <- plotEL(beta1.seq, logel.seq[1,], beta0[1], NA, expression(beta[0]))
logelmode2 <- plotEL(beta2.seq, logel.seq[2,], beta0[2], NA, expression(beta[1]))

# 3-d plot
Beta.seq <- as.matrix(expand.grid(beta1.seq, beta2.seq))
counter <- 0
system.time({
  logel.mat <- apply(Beta.seq, 1, function(bb) {
    counter <<- counter+1
    if (counter %% 100 == 0) message("counter = ", counter)
    -mr.neglogEL_R(y,X,c(bb[1],bb[2]))
  })
})
logel.mat <- matrix(logel.mat, numpoints, numpoints)
logel.mat[is.infinite(logel.mat)] <- NaN
anyNA(logel.mat)
plot3d(seq.x=beta1.seq, seq.y=beta2.seq, seq.z=t(logel.mat))

# mode quadrature approximation for CI
# library(mvtnorm)
# cbind(qmvnorm(p=0.025, mean=nlmout$estimate, sigma=sigma)$quantile,
#       qmvnorm(p=0.975, mean=nlmout$estimate, sigma=sigma)$quantile)
est <- nlmout$estimate
sigma <- solve(nlmout$hessian)
CI1 <- qnorm(p=c(0.025,0.975),mean=nlmout$estimate,sd=sqrt(sigma[1,1]))
CI2 <- qnorm(p=c(0.025,0.975),mean=nlmout$estimate,sd=sqrt(sigma[2,2]))
names(CI1) <- c('L','U')
names(CI2) <- c('L','U')
CI <- rbind(CI1,CI2)
cbind(CI,est)

# try 95% coverage prob here
nrep <- 200
n <- 400
p <- 2
beta_intercept <- 1
beta_slope <- 1.5
beta0 <- c(beta_intercept, beta_slope)
cover <- c(0,0)
len <- c(0,0)
for (ii in 200:nrep) {
  if (ii %% 10 == 0 ) message("ii = ", ii)
  X0 <- matrix(rep(1,n),n,1)
  X1 <- matrix(rnorm(n),n,1)
  X <- cbind(X0,X1)
  eps <- rnorm(n) # N(0,1) error term
  y <- 1 + c(X1 %*% beta_slope) + eps 
  # nlmout <- tryCatch(expr = nlm(mr.neglogEL_R,c(rnorm(1,mean=beta0[1],sd=0.1),
  #                                               rnorm(1,mean=beta0[2],sd=0.1)),y=y,X=X,hessian=TRUE),
  #                    error = function(e) {message("error"); next})
  nlmout <- nlm(mr.neglogEL_R,c(rnorm(1,mean=beta0[1],sd=0.1),
                                rnorm(1,mean=beta0[2],sd=0.1)),y=y,X=X,hessian=TRUE)
  est <- nlmout$estimate
  sigma <- solve(nlmout$hessian)
  CI1 <- qnorm(p=c(0.025,0.975),mean=nlmout$estimate,sd=sqrt(sigma[1,1]))
  CI2 <- qnorm(p=c(0.025,0.975),mean=nlmout$estimate,sd=sqrt(sigma[2,2]))
  if (beta0[1] >= CI1[1] && beta0[1] <= CI1[2]) cover[1] <- cover[1] + 1
  if (beta0[2] >= CI2[1] && beta0[2] <= CI2[2]) cover[2] <- cover[2] + 1
  len[1] <- len[1] + (CI1[2]-CI1[1])
  len[2] <- len[2] + (CI2[2]-CI2[1])
}
cover/nrep
len/nrep

# 3-d problem
n <- 200
p <- 2
X0 <- matrix(rep(1,n),n,1)
X1 <- matrix(rnorm(2*n),n,2)
X <- cbind(X0,X1)
eps <- rnorm(n) # N(0,1) error term
beta_intercept <- 1
beta_slope <- c(1.5,-1.5)
y <- 1 + c(X1 %*% beta_slope) + eps 
beta0 <- c(beta_intercept, beta_slope)

mr.neglogEL_R(y,X,beta0)

nlm(mr.neglogEL_R,beta0*1.05,y=y,X=X)

## qr ##

# 1-d problem 
n <- 200
mu0 <- 1
X <- matrix(rep(1,n), n, 1) # each row of X is one observation
eps <- rnorm(n) # N(0,1) error term
y <- c(X * mu0) + eps
tau <- 0.75

qr.neglogEL.smooth_R(y,X,tau,mu0)

numpoints <- 100
mu.seq <- seq(-.5+mu0+qnorm(tau),.5+mu0+qnorm(tau),length.out = numpoints)
logel.seq <- rep(NA,numpoints)
grad.seq <- rep(NA,numpoints)
for (ii in 1:numpoints) {
  temp <- -qr.neglogEL.smooth_R(y,X,tau,mu.seq[ii],s=10)
  logel.seq[ii] <- temp
  grad.seq[ii] <- attributes(temp)$gradient
}
logelmode <- plotEL(mu.seq, logel.seq, mu0+qnorm(tau), quantile(y,tau), expression(mu))

(nlmout <- nlm(qr.neglogEL.smooth_R,1,y=y,X=X,tau=tau,s=10))
plot(mu.seq,grad.seq,type='l')
abline(h=0,col='blue')
abline(v=nlmout$estimate,col='blue')

# 2-d problem
n <- 500
p <- 2
X0 <- matrix(rep(1,n),n,1)
X1 <- matrix(rnorm(n),n,1)
X <- cbind(X0,X1)
# eps <- rnorm(n) # N(0,1) error term
eps <- gen_eps(n,dist="nct",df=10,ncp=-1)
beta_intercept <- 0.5
beta_slope <- 1
y <- 1 + c(X1 %*% beta_slope) + eps 
beta0 <- c(beta_intercept, beta_slope)
tau <- 0.75
gen_eps(1,dist="nct",df=10,ncp=-1,tau=tau)$nu0+beta0[1]
beta.hat <- coef(rq(y~X-1),tau)
beta.hat
qr.neglogEL_R(y,X,tau,beta0,s = 10)

nlm(qr.neglogEL.smooth_R,beta.hat,y=y,X=X,tau=tau)

# 3-d problem
n <- 200
p <- 2
X0 <- matrix(rep(1,n),n,1)
X1 <- matrix(rnorm(2*n),n,2)
X <- cbind(X0,X1)
eps <- rnorm(n) # N(0,1) error term
beta_intercept <- 1
beta_slope <- c(1.5,-1.5)
y <- 1 + c(X1 %*% beta_slope) + eps 
beta0 <- c(beta_intercept, beta_slope)
tau <- 0.75

qr.neglogEL_R(y,X,tau,beta0)

nlm(qr.neglogEL_R,beta0*1.05,y=y,X=X,tau=tau)

## mr.cens
# TODO: the gradient is not correct right now.

# 1-d problem
n <- 200
mu0 <- 1
X <- matrix(rep(1,n), n, 1) # each row of X is one observation
eps <- rnorm(n) # N(0,1) error term

# random censoring
# yy <- c(X * mu0) + eps
# cc <- rnorm(n,mean=2.5,sd=1)
# deltas <- yy<=cc
# y <- yy
# sum(1-deltas)/n
# y[as.logical(1-deltas)] <- cc[as.logical(1-deltas)]
cc <- rnorm(n,mean=1.35,sd=1)
deltas <- eps<=cc
sum(1-deltas)/n
eps[as.logical(1-deltas)] <- cc[as.logical(1-deltas)]
y <- c(X * mu0) + eps

numpoints <- 100
mu.seq <- seq(-.5+mu0,.5+mu0,length.out = numpoints)
logel.seq <- rep(NA,numpoints)
grad.seq <- rep(NA,numpoints)
for (ii in 1:numpoints) {
  if (ii %% 10 == 0) message("ii = ", ii)
  temp <- -mr_cens.neglogEL_R(y,X,deltas,mu.seq[ii])
  logel.seq[ii] <- temp
  # grad.seq[ii] <- attributes(temp)$gradient
}
logelmode <- plotEL(mu.seq, logel.seq, mu0, mean(y), expression(mu))
# plot(grad.seq,type='l')
# abline(h=0,col='blue')

mr_cens.neglogEL_R(y,X,deltas,1.045277)

nlm(mr_cens.neglogEL_R,mean(y),y=y,X=X,deltas=deltas,hessian=TRUE)

# 2-d problem 
n <- 200
p <- 2
X1 <- matrix(rnorm(n),n,1)
# X1 <- matrix(sample(15:25,n,replace = TRUE),n,1)
X <- cbind(1,X1)
# eps <- rnorm(n) # N(0,1) error term

# dist is one of "norm","t","chisq","lnorm"
eps <- gen_eps(n, dist = "norm", df = NULL)

beta_I <- 0.5
beta_S <- 1
beta0 <- c(beta_I, beta_S)
# plot(X1,y,cex=0.3)

# random censoring
# yy <- beta_I + c(X1 %*% beta_S) + eps 
# cc <- rnorm(n,mean=2,sd=1)
# deltas <- yy<=cc
# y <- yy
# sum(1-deltas)/n
# y[as.logical(1-deltas)] <- cc[as.logical(1-deltas)]

# random censoring
cc <- rnorm(n,mean=1.35,sd=1)
deltas <- eps<=cc
sum(1-deltas)/n
eps[as.logical(1-deltas)] <- cc[as.logical(1-deltas)]
y <- beta_I + c(X1 %*% beta_S) + eps 
plot(X1,y,cex=0.3)

mr_cens.neglogEL_R(y,X,deltas,beta0)

nlm(mr_cens.neglogEL_R,c(0.35,1.2),y=y,X=X,deltas=deltas,hessian = TRUE)

# ---- try nlm ----
negnorm <- function(x) {
  return(-dnorm(x, mean=1))
}
nlm(negnorm,1)

## qrls.smooth

# check if the hessian is close to what nlm return

n <- 200
tau <- 0.75
beta0 <- c(0.5,1)
gamma0 <- -0.5
sig20 <- 1

X <- cbind(rep(1,n),rnorm(n))
Z <- matrix(rnorm(1*n),n,1)
genout <- gen_eps(n, dist = "norm", tau = tau)
eps <- genout$eps
nu0 <- genout$nu0
theta0 <- c(beta0,gamma0,sig20,nu0)
y <- c(X %*% beta0 + sqrt(sig20)*exp(Z %*% gamma0)*eps)
plot(y~X[,2], cex=0.3)

# initial value
hlmout <- hlm.fit(y, X, cbind(1,Z))
hlmout$conv

beta.hlm <- hlmout$beta
gamma.hlm <- hlmout$gamma[2]*0.5
sig2.hlm <- exp(hlmout$gamma[1])
nu.hlm <- quantile((y-X %*% beta.hlm)*exp(-Z %*% gamma.hlm)/sqrt(sig2.hlm),tau)
(theta.hlm <- c(beta.hlm,gamma.hlm,sig2.hlm,nu.hlm))

(nlmout <- nlm(f=qrls.neglogEL.smooth_R,p=theta.hlm,y=y,X=X,Z=Z,tau=tau,s=10,ret_grad=FALSE,
               check.analyticals = TRUE, print.level = 2))
(theta.nlm <- nlmout$estimate)
nlmout$gradient

qrls.neglogEL.smooth_R(y,X,Z,tau,theta.nlm)

# what if fix nu.. does not work either:(
# opt_bool: T if an entry is to be optimized
# theta.fix: fixed params not to be optimized
qrls.neglogEL.smooth_R_fix <- function(y,X,Z,tau,theta.opt,theta.fix,opt_bool,s=10,ret_grad=TRUE) {
  theta.new <- rep(NA,ncol(X)+ncol(Z)+2)
  theta.new[opt_bool] <- theta.opt
  theta.new[!opt_bool] <- theta.fix
  retval <- qrls.neglogEL.smooth_R(y,X,Z,tau,theta.new,s,ret_grad)
  gradval <- attributes(retval)$gradient
  gradvalmn <- gradval[opt_bool]
  attr(retval, "gradient") <- gradvalmn
  retval
}

qrls.neglogEL.smooth_R_fix(y,X,Z,tau,beta.hlm[1],
                           c(beta.hlm[2],gamma.hlm,sig2.hlm,nu.hlm),
                           opt_bool=c(T,F,F,F,F))
(nlmout <- nlm(f=qrls.neglogEL.smooth_R_fix,p=sig2.hlm,
               theta.fix=c(beta.hlm,gamma.hlm,nu.hlm),opt_bool=c(F,F,F,T,F),
               y=y,X=X,Z=Z,tau=tau,ret_grad=TRUE))

# check as one param moves, how does its gradient move

numpoints <- 300
beta0.seq <- seq(-.5+theta.nlm[1], .5+theta.nlm[1], length.out = numpoints)
logel.beta0 <- rep(NA,numpoints)
beta0.grad <- rep(NA,numpoints)
for (ii in 1:numpoints) {
  theta.temp <- theta.nlm
  theta.temp[1] <- beta0.seq[ii]
  retval <- qrls.neglogEL.smooth_R(y,X,Z,tau,theta.temp)
  logel.beta0[ii] <- retval
  beta0.grad[ii] <- attributes(retval)$gradient[1]
}
plot(beta0.seq,logel.beta0,cex=.2)
plot(beta0.seq,beta0.grad,cex=.2)
abline(h=0,col='red')
abline(v=theta.nlm[1],col='red')

beta1.seq <- seq(-.5+theta.nlm[2], .5+theta.nlm[2], length.out = numpoints)
logel.beta1 <- rep(NA,numpoints)
beta1.grad <- rep(NA,numpoints)
for (ii in 1:numpoints) {
  theta.temp <- theta.nlm
  theta.temp[2] <- beta1.seq[ii]
  retval <- qrls.neglogEL.smooth_R(y,X,Z,tau,theta.temp)
  logel.beta1[ii] <- retval
  beta1.grad[ii] <- attributes(retval)$gradient[2]
}
plot(beta1.seq,logel.beta1,cex=.2)
plot(beta1.seq,beta1.grad,cex=.2)
abline(h=0,col='red')
abline(v=theta.nlm[2],col='red')

gamma.seq <- seq(-.5+theta.nlm[3], .5+theta.nlm[3], length.out = numpoints)
logel.gamma <- rep(NA,numpoints)
gamma.grad <- rep(NA,numpoints)
for (ii in 1:numpoints) {
  theta.temp <- theta.nlm
  theta.temp[3] <- gamma.seq[ii]
  retval <- qrls.neglogEL.smooth_R(y,X,Z,tau,theta.temp)
  logel.gamma[ii] <- retval
  gamma.grad[ii] <- attributes(retval)$gradient[3]
}
plot(gamma.seq,logel.gamma,cex=.2)
plot(gamma.seq,gamma.grad,cex=.2)
abline(h=0,col='red')
abline(v=theta.nlm[3],col='red')

sig2.seq <- seq(-.5+theta.nlm[4], .5+theta.nlm[4], length.out = numpoints)
logel.sig2 <- rep(NA,numpoints)
sig2.grad <- rep(NA,numpoints)
for (ii in 1:numpoints) {
  theta.temp <- theta.nlm
  theta.temp[4] <- sig2.seq[ii]
  retval <- qrls.neglogEL.smooth_R(y,X,Z,tau,theta.temp)
  logel.sig2[ii] <- retval
  sig2.grad[ii] <- attributes(retval)$gradient[4]
}
plot(sig2.seq,logel.sig2,cex=.2)
plot(sig2.seq,sig2.grad,cex=.2)
abline(h=0,col='red')
abline(v=theta.nlm[4],col='red')

nu.seq <- seq(-.5+theta.nlm[5], .5+theta.nlm[5], length.out = numpoints)
logel.nu <- rep(NA,numpoints)
nu.grad <- rep(NA,numpoints)
for (ii in 1:numpoints) {
  theta.temp <- theta.nlm
  theta.temp[5] <- nu.seq[ii]
  retval <- qrls.neglogEL.smooth_R(y,X,Z,tau,theta.temp)
  logel.nu[ii] <- retval
  nu.grad[ii] <- attributes(retval)$gradient[5]
}
plot(nu.seq,logel.nu,cex=.2)
plot(nu.seq,nu.grad,cex=.2)
abline(h=0,col='red')
abline(v=theta.nlm[5],col='red')

