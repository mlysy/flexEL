## 3-d plot of logEL surface 
library(bayesEL)
source("../testthat/el-utils.R")
source("../testthat/el-rfuns.R")
source("../testthat/el-model.R")
source("gen_eps.R")
source("mode-functions.R")
source("smoothEL.R")
library(plotly)

# check how the coordinate works with plot_ly
# plot_ly(x=c(-1,0,1),y=c(2,3,4),z=matrix(c(-4,-3,-2,-1,0,1,2,3,4),3,3), type="surface")

# ---- HMC paper Example 1 ----
hmc.evalG <- function(y, v, theta) {
  G <- matrix(NA,nrow=length(y),ncol=2)
  G[,1] <- y-(theta[1]+theta[2]*v)
  G[,2] <- v*G[,1]
  return(G)
}

hmc.logEL <- function(y,v,theta) {
  G <- hmc.evalG(y,v,theta)
  omegas <- omega.hat(G)
  # indep normal prior N(0,k^2)
  logel <- sum(log(omegas)) + sum(dnorm(theta,mean=0,sd=100,log = TRUE))
  return(logel)
}

n <- 30
theta <- c(0.5,1)
v <- rnorm(n)
e <- rnorm(n)
y <- theta[1] + theta[2]*v + e

numpoints <- 100
theta1.seq <- seq(theta[1]-1.5,theta[1]+1.5,length.out = numpoints)
theta2.seq <- seq(theta[2]-1.5,theta[2]+1.5,length.out = numpoints)
theta.seq <- cbind(theta1.seq,theta2.seq)
logel.seq <- matrix(rep(NA,2*numpoints),2,numpoints)

Theta.seq <- as.matrix(expand.grid(theta1.seq, theta2.seq))
logel.mat <- apply(Theta.seq, 1, function(tt) {
  hmc.logEL(y,v,tt)
})
logel.mat <- matrix(logel.mat, numpoints, numpoints)
# el.mat <- exp(logel.mat - max(logel.mat))

# mode point
modept <- which(logel.mat == max(logel.mat,na.rm = TRUE), arr.ind = TRUE)
c(theta1.seq[modept[1]],theta2.seq[modept[2]])

# linear regression 
lmcoef <- coef(lm(y ~ v))

# 3d surface plot
plot3d(seq.x=theta1.seq, seq.y=theta2.seq, seq.z=t(logel.mat))

# ---- mr ----
n <- 50
p <- 2
X0 <- matrix(rep(1,n),n,1)
X1 <- matrix(rnorm(n),n,1)
X <- cbind(X0,X1)
eps <- rnorm(n) # N(0,1) error term
beta_intercept <- 0.5
beta_slope <- 1.0
y <- 1 + c(X1 %*% beta_slope) + eps 
beta0 <- c(beta_intercept, beta_slope)
plot(X1,y,cex=0.3)

numpoints <- 100
beta1.seq <- seq(beta0[1]-1.5,beta0[1]+1.5,length.out = numpoints)
beta2.seq <- seq(beta0[2]-1.5,beta0[2]+1.5,length.out = numpoints)
logel.seq <- matrix(rep(NA,2*numpoints),2,numpoints)

Beta.seq <- as.matrix(expand.grid(beta1.seq, beta2.seq))
logel.mat <- apply(Beta.seq, 1, function(bb) {
  # attributes(-mr.neglogEL_R(y,X,bb))$gradient[2]
  G <- mr.evalG(y,X,c(bb[1],bb[2]))
  omegas <- omega.hat(G)
  logEL(omegas)
})
logel.mat <- matrix(logel.mat, numpoints, numpoints)
logel.mat[is.infinite(logel.mat)] <- NaN
anyNA(logel.mat)
# el.mat <- exp(logel.mat - max(logel.mat))

# linear regression 
lmcoef <- coef(lm(y ~ X1))

# plot 3d surface
G <- mr.evalG(y,X,beta0)
omegas <- omega.hat(G)
logel.true <- logEL(omegas)
G <- mr.evalG(y,X,lmcoef)
omegas <- omega.hat(G)
logel.est <- logEL(omegas)
plot3d(seq.x=beta1.seq, seq.y=beta2.seq, seq.z=t(logel.mat),
       m1=beta0, logel.m1=logel.true,
       m2=lmcoef, logel.m2=logel.est)
# plot3d(seq.x=beta1.seq, seq.y=beta2.seq, seq.z=t(logel.mat))

# ---- mr.cens (smoothed) ----
n <- 100
p <- 2
X1 <- matrix(rnorm(n),n,1)
# X1 <- matrix(sample(15:25,n,replace = TRUE),n,1)
X <- cbind(1,X1)
# eps <- rnorm(n) # N(0,1) error term

# dist is one of "norm","t","chisq","lnorm"
eps <- gen_eps(n, dist = "norm", df = NULL)

beta_I <- 0.5
beta_S <- 1
# beta_I <- rnorm(1)
# beta_S <- rnorm(1)
yy <- beta_I + c(X1 %*% beta_S) + eps 
beta0 <- c(beta_I, beta_S)
# plot(X1,y,cex=0.3)

# random censoring
cc <- rnorm(n,mean=2,sd=1)
deltas <- yy<=cc
y <- yy
sum(1-deltas)/n
y[as.logical(1-deltas)] <- cc[as.logical(1-deltas)]

numpoints <- 40
beta1.seq <- seq(beta0[1]-1.5,beta0[1]+1.5,length.out = numpoints)
beta2.seq <- seq(beta0[2]-1.5,beta0[2]+1.5,length.out = numpoints)
beta.seq <- cbind(beta1.seq,beta2.seq)
logel.seq <- matrix(rep(NA,2*numpoints),2,numpoints)

# TODO: try using the accelerated logEL with adjusted G
Beta.seq <- as.matrix(expand.grid(beta1.seq, beta2.seq))
adjust <- FALSE
logel.mat <- apply(Beta.seq, 1, function(bb) {
  G <- mr.evalG(y,X,c(bb[1],bb[2]))
  epsilons <- y - c(X %*% c(bb[1],bb[2]))
  if (adjust) {
    G <- adjG_R(G)
    omegas <- omega.hat_R(G,deltas,epsilons,adjust)$omegas
    logEL(omegas,epsilons,deltas,adjust)
  }
  else {
    omegas <- omega.hat.EM.smooth_R(G,deltas,epsilons)$omegas # smoothed version
    logEL(omegas,epsilons,deltas)
  }
})logel.mat <- matrix(logel.mat, numpoints, numpoints)
logel.mat[is.infinite(logel.mat)] <- NaN
anyNA(logel.mat)
# el.mat <- exp(logel.mat - max(logel.mat))

# linear regression 
lmcoef <- coef(lm(y ~ X1))

# plot 3d surface
G <- mr.evalG(y,X,beta0)
epsilons <- y - c(X %*% beta0)
omegas <- omega.hat(G,deltas,epsilons)
logel.true <- logEL(omegas,epsilons,deltas)
G <- mr.evalG(y,X,lmcoef)
epsilons <- y - c(X %*% lmcoef)
omegas <- omega.hat(G,deltas,epsilons)
logel.est <- logEL(omegas,epsilons,deltas)
plot3d(seq.x=beta1.seq, seq.y=beta2.seq, seq.z=t(logel.mat),
       m1=beta0, logel.m1=logel.true,
       m2=lmcoef, logel.m2=logel.est)

# ---- qr ----
n <- 100
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

Beta.seq <- as.matrix(expand.grid(beta1.seq, beta2.seq))
logel.mat <- apply(Beta.seq, 1, function(bb) {
  G <- qr.evalG(y,X,alpha,matrix(c(bb[1],bb[2]),2,1))
  omegas <- omega.hat(G)
  logEL(omegas)
})
logel.mat <- matrix(logel.mat, numpoints, numpoints)
logel.mat[is.infinite(logel.mat)] <- NaN
anyNA(logel.mat)
# el.mat <- exp(logel.mat - max(logel.mat))
# logel.marg <- log(cbind(beta1 = rowSums(el.mat), beta2 = colSums(el.mat)))

# quantile regression
library(quantreg)
qrcoef <- coef(rq(y~X1,tau=alpha))

# plot 3d surface
G <- qr.evalG(y,X,alpha,beta0)
omegas <- omega.hat(G)
logel.true <- logEL(omegas)
G <- qr.evalG(y,X,alpha,qrcoef)
omegas <- omega.hat(G)
logel.est <- logEL(omegas)
plot3d(seq.x=theta1.seq, seq.y=theta2.seq, seq.z=t(logel.mat),
       m1=beta0, logel.m1=logel.true,
       m2=qrcoef, logel.m2=logel.est)

# ---- qr.cens ----
n <- 100
p <- 2
alpha <- 0.5
X0 <- matrix(rep(1,n),n,1)
# X1 <- matrix(seq(-2,2,length.out = n),n,1)
X1 <- matrix(rnorm(n),n,1)
X <- cbind(X0,X1)

# eps <- rnorm(n) # N(0,1) error term
genout <- gen_eps(n, dist = "norm", df = NULL, tau = alpha)
eps <- genout$eps
nu0 <- genout$nu0

beta_I <- 0.5
beta_S <- 1
yy <- c(X1 %*% beta_S) + eps 
plot(X1,yy,cex=0.3)
beta0 <- c(beta_I, beta_S)
beta0

# random censoring
cc <- rnorm(n,mean=1.5,sd=1)
deltas <- yy<=cc
y <- yy
sum(1-deltas)/n
y[as.logical(1-deltas)] <- cc[as.logical(1-deltas)]

numpoints <- 100
beta1.seq <- seq(beta0[1]-1.5,beta0[1]+1.5,length.out = numpoints)
beta2.seq <- seq(beta0[2]-1.5,beta0[2]+1.5,length.out = numpoints)
# theta.seq <- cbind(theta1.seq,theta2.seq)
logel.seq <- matrix(rep(NA,2*numpoints),2,numpoints)

Beta.seq <- as.matrix(expand.grid(beta1.seq, beta2.seq))
logel.mat <- apply(Beta.seq, 1, function(bb) {
  G <- qr.evalG(y,X,alpha,matrix(c(bb[1],bb[2]),2,1))
  epsilons <- y - c(X %*% matrix(c(bb[1],bb[2]),2,1))
  omegas <- omega.hat(G,deltas,epsilons)
  logEL(omegas,epsilons,deltas)
})
logel.mat <- matrix(logel.mat, numpoints, numpoints)
logel.mat[is.infinite(logel.mat)] <- NaN
anyNA(logel.mat)
# el.mat <- exp(logel.mat - max(logel.mat))
# logel.marg <- log(cbind(beta1 = rowSums(el.mat), beta2 = colSums(el.mat)))

# quantile regression
library(quantreg)
qrcoef <- coef(rq(y~X1,tau=alpha))

# plot 3d surface
G <- qr.evalG(y,X,alpha,beta0)
epsilons <- y - c(X %*% beta0)
omegas <- omega.hat(G,deltas,epsilons)
logel.true <- logEL(omegas,epsilons,deltas)
G <- qr.evalG(y,X,alpha,qrcoef)
epsilons <- y - c(X %*% qrcoef)
omegas <- omega.hat(G,deltas,epsilons)
logel.est <- logEL(omegas,epsilons,deltas)
plot3d(seq.x=theta1.seq, seq.y=theta2.seq, seq.z=t(logel.mat),
       m1=beta0, logel.m1=logel.true,
       m2=qrcoef, logel.m2=logel.est)

# ---- smoothed qr ----
# 2-d problem
n <- 200
p <- 2
X0 <- matrix(rep(1,n),n,1)
X1 <- matrix(rnorm(n),n,1)
X <- cbind(X0,X1)
tau <- 0.75
genout <- gen_eps(n, dist = "norm", df = NULL, tau = tau)
eps <- genout$eps
nu0 <- genout$nu0
beta_intercept <- 0.5
beta_slope <- 1
y <- 1 + c(X1 %*% beta_slope) + eps 
beta0 <- c(beta_intercept+nu0, beta_slope)

numpoints <- 100
beta1.seq <- seq(beta0[1]-1.5,beta0[1]+1.5,length.out = numpoints)
beta2.seq <- seq(beta0[2]-1.5,beta0[2]+1.5,length.out = numpoints)
beta.seq <- cbind(beta1.seq,beta2.seq)
logel.seq <- matrix(rep(NA,2*numpoints),2,numpoints)

Beta.seq <- as.matrix(expand.grid(beta1.seq, beta2.seq))
logel.mat <- apply(Beta.seq, 1, function(bb) {
  G <- qr.evalG.smooth_R(y,X,tau,bb,s=10)
  omegas <- omega.hat(G)
  logEL(omegas)
})
logel.mat <- matrix(logel.mat, numpoints, numpoints)
logel.mat[is.infinite(logel.mat)] <- NaN
anyNA(logel.mat)

# quantile regression
library(quantreg)
qrcoef <- coef(rq(y~X1,tau=alpha))

# plot 3d surface
G <- qr.evalG.smooth_R(y,X,tau,beta0,s=10)
omegas <- omega.hat(G)
logel.true <- logEL(omegas)
G <- qr.evalG.smooth_R(y,X,tau,qrcoef,s=10)
omegas <- omega.hat(G)
logel.est <- logEL(omegas)
plot3d(seq.x=beta1.seq, seq.y=beta2.seq, seq.z=t(logel.mat),
       m1=beta0, logel.m1=logel.true,
       m2=qrcoef, logel.m2=logel.est)

# ---- EL + censoring: how does the support change when % changes ----
# mr.cens 
n <- 100
p <- 2
X1 <- matrix(rnorm(n),n,1)
# X1 <- matrix(sample(15:25,n,replace = TRUE),n,1)
X <- cbind(1,X1)
# eps <- rnorm(n) # N(0,1) error term

# dist is one of "norm","t","chisq","lnorm"
eps <- gen_eps(n, dist = "norm", df = NULL)

beta_I <- 0.5
beta_S <- 1
# beta_I <- rnorm(1)
# beta_S <- rnorm(1)
yy <- beta_I + c(X1 %*% beta_S) + eps 
beta0 <- c(beta_I, beta_S)
# plot(X1,y,cex=0.3)

# random censoring
cc <- rnorm(n,mean=0,sd=1)
deltas <- yy<=cc
sum(1-deltas)/n
deltas.save <- deltas

# select a subset of cc so that the censored values keep the same
deltas <- deltas.save
deltas[deltas==0][41:50] <- TRUE

y <- yy
y[as.logical(1-deltas)] <- cc[as.logical(1-deltas)]

numpoints <- 100
beta1.seq <- seq(beta0[1]-1.5,beta0[1]+1.5,length.out = numpoints)
beta2.seq <- seq(beta0[2]-1.5,beta0[2]+1.5,length.out = numpoints)
beta.seq <- cbind(beta1.seq,beta2.seq)
logel.seq <- matrix(rep(NA,2*numpoints),2,numpoints)

Beta.seq <- as.matrix(expand.grid(beta1.seq, beta2.seq))
adjust <- FALSE
logel.mat <- apply(Beta.seq, 1, function(bb) {
  G <- mr.evalG(y,X,c(bb[1],bb[2]))
  if (adjust) G <- adjG_R(G)
  epsilons <- y - c(X %*% c(bb[1],bb[2]))
  if (adjust) {
    omegas <- omega.hat_R(G,deltas,epsilons,adjust)
    logEL_R(omegas,epsilons,deltas,adjust)
  }
  else {
    omegas <- omega.hat(G,deltas,epsilons)
    logEL(omegas,epsilons,deltas)
  }
})
logel.mat <- matrix(logel.mat, numpoints, numpoints)
logel.mat[is.infinite(logel.mat)] <- NaN
anyNA(logel.mat)
# el.mat <- exp(logel.mat - max(logel.mat))

# linear regression 
lmcoef <- coef(lm(y ~ X1))

# plot 3d surface
G <- mr.evalG(y,X,beta0)
epsilons <- y - c(X %*% beta0)
omegas <- omega.hat(G,deltas,epsilons)
logel.true <- logEL(omegas,epsilons,deltas)
G <- mr.evalG(y,X,lmcoef)
epsilons <- y - c(X %*% lmcoef)
omegas <- omega.hat(G,deltas,epsilons)
logel.est <- logEL(omegas,epsilons,deltas)
plot3d(type="surface", seq.x=beta1.seq, seq.y=beta2.seq, seq.z=t(logel.mat),
       m1=beta0, logel.m1=logel.true,
       m2=lmcoef, logel.m2=logel.est)
plot3d(type="contour", seq.x=beta1.seq, seq.y=beta2.seq, seq.z=t(logel.mat),
       title="mr, n=100, 50% censored")


