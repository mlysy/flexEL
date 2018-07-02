# source("../../R/hlm.R")
library(bayesEL)
library(numDeriv)
source("gen_eps.R")
source("../testthat/el-utils.R")
source("../testthat/el-rfuns.R")
source("../testthat/el-model.R")
source("mode-functions.R")
source("em_acceleration.R")

# ---- mr: modeEL v.s. bayesEL ----
# dimension
n <- 200
p <- 2
X1 <- rnorm(n)
X <- cbind(1,X1)

# parameters
beta0 <- c(1,0.5)
beta0

# dist is one of "norm","t","chisq","lnorm"
eps <- gen_eps(n, dist = "chisq", df = 5)

# response
y <- c(X %*% beta0) + eps
plot(y~X1, cex=.3)

# fit lm
lmcoeff <- coef(lm(y~X1))
lmcoeff
abline(a=lmcoeff[1],b=lmcoeff[2])

# by mode finding 
nlmout <- nlm(mr.neglogEL_R,lmcoeff,y=y,X=X,hessian=TRUE)
nlmout$estimate

# by full Bayesian
nsample <- 10000
nburn <- 2000
belout <- mr.post_adapt(y,X,nsample,nburn,lmcoeff,c(0.1,0.1))
cbind(mean(belout$beta_chain[1,]),mean(belout$beta_chain[2,]))

# compare MSE of estimates: mode finding is better 
# but it seems alwyas coincide with lm result
nrep <- 200
mod.ms <- c(0,0)
bel.ms <- c(0,0)
for (ii in 1:nrep) {
  if (ii %% 10 == 0) message("ii = ", ii)
  X1 <- rnorm(n)
  X <- cbind(1,X1)
  eps <- gen_eps(n, dist = "lnorm", df = NULL)
  y <- c(X %*% beta0) + eps
  # fit lm
  lmcoeff <- coef(lm(y~X1))
  # mode finding
  nlmout <- nlm(mr.neglogEL_R,lmcoeff,y=y,X=X,hessian=TRUE)
  mod.ms <- mod.ms + (nlmout$estimate-beta0)^2
  # full bayesian
  belout <- mr.post_adapt(y,X,nsample,nburn,lmcoeff,c(0.1,0.1))
  bel.ms <- bel.ms + (c(mean(belout$beta_chain[1,]),mean(belout$beta_chain[2,]))-beta0)^2
}
mod.ms/nrep # n=100: 0.009226374 0.009351438; n=200: 0.005084070 0.005386566
bel.ms/nrep # n=100: 0.01780397 0.01808397; n=200: 0.008233900 0.008418917

# ---- mrls ----

# check MSE 

# ---- qr cens ls ---- 

# dimension
n <- 200
p <- 2
q <- 1

# covaraites
# X1 <- rnorm(n)
X1 <- sample(-2:2,n,replace = TRUE)
# X2 <- rbinom(n,1,0.3)
# X <- cbind(X1,X2)
X <- cbind(1,X1)
Z <- matrix(X1)

# parameters
beta0 <- c(1,0.5)
beta0
gamma0 <- -.5
gamma0
sig20 <- 1
sig20

# quantile level
tau <- 0.75

# dist is one of "norm","t","chisq","lnorm"
genout <- gen_eps(n, dist = "norm", df = NULL, tau = tau)
eps <- genout$eps
nu0 <- genout$nu0

# responses and censorings
yy <- c(X %*% beta0) + sqrt(sig20)*exp(0.5*X1*gamma0)*eps
cc <- rnorm(n, mean=3)
delta <- yy <= cc
sum(1-delta)/length(delta)
y <- yy
y[!delta] <- cc[!delta]

# ---- plot the lines ----
plot(y=y,x=X1,cex=.3)
abline(a=beta0[1],b=beta0[2],col='red')
abline(a=mean(theta_chain[1,]),b=mean(theta_chain[2,]),col='blue')
abline(a=beta.hlm[1],b=beta.hlm[2],col='green')

# plot conditional likelihood
numpoints <- 100
deltas <- delta
adjust <- FALSE

beta.seq1 <- seq(-.5+beta0[1],.5+beta0[1],length.out = numpoints)
logel.seq1 <- rep(NA,numpoints)
for (ii in 1:numpoints) {
  message("ii = ", ii)
  G <- qrls.evalG(y,X,Z,tau,c(beta.seq1[ii],beta0[2]),0.5*gamma0,sig20,nu0)
  if (adjust) G <- adjG_R(G)
  epsilons <- evalEpsilonsLS(y,X,Z,c(beta.seq1[ii],beta0[2]),0.5*gamma0,sig20)
  if (adjust) {
    omegas <- omega.hat_R(G,deltas,epsilons,adjust)
    logel.seq1[ii] <- logEL_R(omegas,epsilons,deltas,adjust)
  }
  else {
    omegas <- omega.hat(G,deltas,epsilons)
    logel.seq1[ii] <- logEL(omegas,epsilons,deltas)
  }
}
logelmode1 <- plotEL(beta.seq1, logel.seq1, beta0[1], NA, expression(beta[0]))

beta.seq2 <- seq(-.5+beta0[2],.5+beta0[2],length.out = numpoints)
logel.seq2 <- rep(NA,numpoints)
for (ii in 1:numpoints) {
  message("ii = ", ii)
  G <- qrls.evalG(y,X,Z,tau,c(beta0[1],beta.seq2[ii]),0.5*gamma0,sig20,nu0)
  if (adjust) G <- adjG_R(G)
  epsilons <- evalEpsilonsLS(y,X,Z,c(beta0[1],beta.seq2[ii]),0.5*gamma0,sig20)
  if (adjust) {
    omegas <- omega.hat_R(G,deltas,epsilons,adjust)
    logel.seq2[ii] <- logEL_R(omegas,epsilons,deltas,adjust)
  }
  else {
    omegas <- omega.hat(G,deltas,epsilons)
    logel.seq2[ii] <- logEL(omegas,epsilons,deltas)
  }
}
logelmode2 <- plotEL(beta.seq2, logel.seq2, beta0[2], NA, expression(beta[1]))

gamma.seq <- seq(-.5+gamma0/2,.5+gamma0/2,length.out = numpoints)
logel.seq3 <- rep(NA,numpoints)
for (ii in 1:numpoints) {
  message("ii = ", ii)
  G <- qrls.evalG(y,X,Z,tau,beta0,gamma.seq[ii],sig20,nu0)
  epsilons <- evalEpsilonsLS(y,X,Z,beta0,gamma.seq[ii],sig20)
  omegas <- omega.hat(G,deltas,epsilons)
  logel.seq3[ii] <- logEL(omegas,epsilons,deltas)
}
logelmode3 <- plotEL(gamma.seq, logel.seq3, 0.5*gamma0, NA, expression(gamma[1]))

sig2.seq <- seq(-.5+sig20,.5+sig20,length.out = numpoints)
logel.seq4 <- rep(NA,numpoints)
for (ii in 1:numpoints) {
  message("ii = ", ii)
  G <- qrls.evalG(y,X,Z,tau,beta0,0.5*gamma0,sig2.seq[ii],nu0)
  epsilons <- evalEpsilonsLS(y,X,Z,beta0,0.5*gamma0,sig2.seq[ii])
  omegas <- omega.hat(G,deltas,epsilons)
  logel.seq4[ii] <- logEL(omegas,epsilons,deltas)
}
logelmode4 <- plotEL(sig2.seq, logel.seq4, sig20, NA, expression(sigma^2))

nu.seq <- seq(-.5+nu0,.5+nu0,length.out = numpoints)
logel.seq5 <- rep(NA,numpoints)
for (ii in 1:numpoints) {
  message("ii = ", ii)
  G <- qrls.evalG(y,X,Z,tau,beta0,0.5*gamma0,sig20,nu.seq[ii])
  epsilons <- evalEpsilonsLS(y,X,Z,beta0,0.5*gamma0,sig20)
  omegas <- omega.hat(G,deltas,epsilons)
  logel.seq5[ii] <- logEL(omegas,epsilons,deltas)
}
logelmode5 <- plotEL(nu.seq, logel.seq5, nu0, NA, expression(nu))

# ---- use hlm ----
hlmout <- hlm(y,delta,X,cbind(1,Z),rel_tol = 1e-5)
hlmout$conv
beta.hlm <- hlmout$coef$beta
gamma.hlm <- 0.5*hlmout$coef$gamma[2]
sig2.hlm <- exp(hlmout$coef$gamma[1])
nu.hlm <- quantile((y-X %*% beta.hlm)*exp(-Z %*% gamma.hlm)/sqrt(sig2.hlm),tau)
# quick check 
rbind(true = beta0, est = beta.hlm) # beta
rbind(true = gamma0/2, est = gamma.hlm) # gamma
rbind(true = sig20, est = sig2.hlm) # beta
rbind(true = nu0, est = nu.hlm) # beta

theta.boot <- qrls_cens.boot_R(y,X,Z,deltas,
                               tau,beta.hlm,gamma.hlm,sig2.hlm,nu.hlm,
                               nboot=1000)
# confidence intervals
CI.beta <- sapply(1:length(beta.hlm), function(ii) {
  qrls_cens.bootCI_R(theta.hat = beta.hlm[ii], theta.boot = theta.boot$beta.boot[,ii])
})
CI.beta
CI.gamma <- sapply(1:length(gamma.hlm), function(ii) {
  qrls_cens.bootCI_R(theta.hat = gamma.hlm[ii], theta.boot = theta.boot$gamma.boot[,ii])
})
CI.gamma
CI.sig2 <- qrls_cens.bootCI_R(theta.hat = sig2.hlm, theta.boot = theta.boot$sig2.boot)
CI.sig2
CI.nu <- qrls_cens.bootCI_R(theta.hat = nu.hlm, theta.boot = theta.boot$nu.boot)
CI.nu

# theta.ci <- qrls_cens.bootCI_R(list(beta.hat=beta.hlm,
#                                     gamma.hat=gamma.hlm,
#                                     sig2.hat=sig2.hlm,
#                                     nu.hat=nu.hlm),
#                                theta.boot)

# ---- use bayesEL ----
nsamples <- 5000
nburn <- 1000
betaInit <- beta.hlm
gammaInit <- gamma.hlm
sig2Init <- sig2.hlm
nuInit <- nu.hlm
mwgSd <- rep(0.1,5)
rvDoMcmc <- rep(1,5)
system.time(
  belout <- qrls_cens.post_adapt(y,X,Z,delta,tau,nsamples,nburn,
                                 betaInit,gammaInit,sig2Init,nuInit,
                                 mwgSd,rvDoMcmc)
)
theta_chain <- belout$theta_chain
theta_accept <- belout$paccept
theta_accept

# mrlsout <- postCens_R(Gfun=mrls.evalG_R,nThe=4,nBet=2,nGam=1,
#                       y=y,X=X,Z=matrix(X1),deltas=delta,
#                       thetaInit=c(betaInit,gammaInit,sig2Init),
#                       nsamples=nsamples,nburn=nburn,
#                       mwgSds=mwgSd,adjust = TRUE)
# theta_chain <- mrlsout$theta_chain
# theta_accept <- mrlsout$paccept
# theta_accept

# mixing of the chains
plot(theta_chain[5,],type = 'l')

hist(theta_chain[1,],breaks=50,freq=FALSE,
     xlab = expression(beta[1]), main='')
abline(v=beta0[1],col='red')
abline(v=mean(theta_chain[1,]),col='blue')
abline(v=beta.hlm,col='green')
legend('topright',legend=c(expression('true param'),
                           expression('sample mean'),
                           expression('hlm sol')),
       lty = c(1,1), col = c('red','blue','green'), cex = 0.6)

hist(theta_chain[2,],breaks=50,freq=FALSE,
     xlab = expression(beta[2]), main='')
abline(v=beta0[2],col='red')
abline(v=mean(theta_chain[2,]),col='blue')
abline(v=hlmout$coef$beta[2],col='green')
legend('topright',legend=c(expression('true param'),
                           expression('sample mean'),
                           expression('hlm sol')),
       lty = c(1,1), col = c('red','blue','green'), cex = 0.6)

hist(theta_chain[3,],breaks=50,freq=FALSE,
     xlab = expression(gamma[1]), main='')
abline(v=gamma0[1]/2,col='red')
abline(v=mean(theta_chain[3,]),col='blue')
abline(v=hlmout$coef$gamma[2]/2,col='green')
legend('topright',legend=c(expression('true param'),
                           expression('sample mean'),
                           expression('hlm sol')),
       lty = c(1,1), col = c('red','blue','green'), cex = 0.6)

hist(theta_chain[4,],breaks=50,freq=FALSE,
     xlab = expression(gamma[2]), main='')
abline(v=sig20,col='red')
abline(v=mean(theta_chain[4,]),col='blue')
abline(v=sig2.hlm,col='green')
legend('topright',legend=c(expression('true param'),
                           expression('sample mean'),
                           expression('hlm sol')),
       lty = c(1,1), col = c('red','blue','green'), cex = 0.6)

hist(theta_chain[5,],breaks=50,freq=FALSE,
     xlab = expression(sigma^2), main='')
abline(v=nu0,col='red')
abline(v=mean(theta_chain[5,]),col='blue')
abline(v=nu.hlm,col='green')
legend('topright',legend=c(expression('true param'),
                           expression('sample mean'),
                           expression('hlm sol')),
       lty = c(1,1), col = c('red','blue','green'), cex = 0.6)

# --- test the accelerated EM ---

# with random data
n <- 10
p <- 2
max_iter <- 100
rel_tol <- 1e-5
G <- matrix(rnorm(n*p), n, p)
deltas <- rep(1,n)
numcens <- sample(round(n/2),1)
censinds <- sample(n,numcens)
deltas[censinds] <- 0
epsilons <- rnorm(n)
# omegahat <- omega.hat.EM_Acc_R(G, deltas, epsilons, adjust = FALSE, 
#                         max_iter = max_iter, rel_tol = rel_tol, verbose = TRUE)
omegahat <- omega.hat.EM_R(G, deltas, epsilons, adjust = FALSE, 
                           max_iter = max_iter, rel_tol = rel_tol, verbose = TRUE, dbg= TRUE)

# fit the EM sequence with a curve dbg=TRUE above
idx <- 1/(1:length(omegahat$logels))
pred <- predict(lm(omegahat$logels ~ idx + I(idx^2)))
plot(pred, ylim=range(omegahat$logels,pred))
points(omegahat$logels, col='red')

# fit first 5 points and predict the rest
nfit <- 9
idx <- 1/(1:length(omegahat$logels))
idx_fit <- idx[1:nfit]
idx_pre <- idx[(nfit+1):length(idx)]
fit <- lm(omegahat$logels[1:nfit] ~ idx_fit + I(idx_fit^2))
pred <- predict(fit, newdata = data.frame(idx_fit=idx))
plot(pred, ylim=range(omegahat$logels,pred))
points(omegahat$logels, col='red')

# check EM convergence rate
nrep <- 100
n <- 100
p <- 2
max_iter <- 100
rel_tol <- 1e-5
niter <- rep(NA,nrep)
cens <- rep(NA,nrep)
count <- 0
nnotconv <- 0
while (count < nrep) {
  count <- count+1
  if (count %% 10 == 0) message("count = ", count)
  G <- matrix(rnorm(n*p), n, p)
  deltas <- rep(1,n)
  numcens <- sample(n,1)
  if (numcens/n < 0.89 || numcens/n > 0.91) {
    count <- count-1
    next
  }
  else message("censored pct: ", numcens/n)
  censinds <- sample(n,numcens)
  deltas[censinds] <- 0
  epsilons <- rnorm(n)
  omegahat <- omega.hat.EM_R(G, deltas, epsilons, adjust = FALSE, 
                             max_iter = max_iter, rel_tol = rel_tol, verbose = FALSE, dbg = TRUE)
  # message("omegahat$nIter = ", omegahat$nIter)
  if (!omegahat$conv) {
    nnotconv <- nnotconv + 1
    count <- count-1
    next
  }
  cens[count] <- numcens/n
  niter[count] <- omegahat$nIter
}
mean(cens) # 0.1006, 0.1995, 0.2997, 0.4001, 0.4986, 0.6, 0.7, 0.7999, 0.8982
mean(niter) # 5.75, 8.94, 12.91, 17.24, 23.73, 32.78, 44.55, 83.57
nnotconv
# plot them
censvec <- c(0, 0.1006, 0.1995, 0.2997, 0.4001, 0.4986, 0.6, 0.7, 0.7999, 0.8982)
nitervec <- c(1, 5.75, 8.94, 12.91, 17.24, 23.73, 32.78, 44.55, 64.2, 83.57)
plot(censvec,nitervec, xlab="censored %", ylab="nIter", cex=.5, 
     main = "rel_tol=1e-5, n=100, p=2, rand G; nrep=100")

# with qrls data
# dimensions
n <- 200 # number of observations
p <- 1
q <- 1

# covariates
X <- matrix(rep(1,n),n,p)
# X <- matrix(rnorm(n),n,p)
# Z <- matrix(rep(1,n),n,q) # Z should not contain an intercept
Z <- matrix(rnorm(n),n,q)

# parameters
beta0 <- rnorm(p)
beta0
gamma0 <- rnorm(q)
gamma0
sig20 <- abs(rnorm(1,mean=1)) # NEW: scale param
sig20

# quantile level
alpha <- 0.75

# dist is one of "norm","t","chisq","lnorm"
genout <- gen_eps(n, dist = "norm", df = NULL, tau = alpha)
eps <- genout$eps
nu0 <- genout$nu0

# resposes
yy <- c(X %*% beta0 + sqrt(sig20)*exp(0.5 * Z %*% gamma0)*eps)
plot(yy,cex=0.3)
# random censoring
cc <- rnorm(n,mean=1,sd=1)
deltas <- yy<=cc
y <- yy
sum(1-deltas)/n
y[as.logical(1-deltas)] <- cc[as.logical(1-deltas)]

epsilons <- evalEpsilonsLS_R(y,X,Z,beta0,gamma0,sig20)
G <- qrls.evalG(y,X,Z,alpha,beta0,gamma0,sig20,nu0)

max_iter <- 100
rel_tol <- 1e-5
# omegahat <- omega.hat.EM_Acc_R(G, deltas, epsilons, adjust = FALSE, 
#                                max_iter = max_iter, rel_tol = rel_tol, verbose = TRUE)
# compare the normal way and the acc way
system.time(
  print(logEL_EMAC_R(G,epsilons,deltas))
)

system.time({
  omegahat <- omega.hat.EM_R(G, deltas, epsilons, adjust = FALSE, 
                             max_iter = max_iter, rel_tol = rel_tol, verbose = TRUE, dbg = TRUE)
  print(logEL_R(omegahat$omegas,epsilons,deltas))
})

# ---- compare the coverage prob of the two methods ----
# repeat the experiment and calculate coverage probilities
con <- file("covprob.log")
sink(con, append=TRUE)
sink(con, append=TRUE, type="message")

ns <- c(100,200,300)
alpha <- 0.75
beta0 <- c(1,0.5)
gamma0 <- -.5
sig20 <- 1

nsamples <- 10000
nburn <- 5000
mwgSds <- rep(0.1,5)
reptimes <- 200
coverProbs.el <- matrix(rep(0,5*length(ns)),5,length(ns))
coverProbs.hlm <- matrix(rep(0,5*length(ns)),5,length(ns))
lengthCI.el <- matrix(rep(0,5*length(ns)),5,length(ns))
lengthCI.hlm <- matrix(rep(0,5*length(ns)),5,length(ns))

# change this according to the eps in the loop
nu0 <- gen_eps(100, dist = "chisq", df = 5, tau = alpha)$nu0
theta0 <- c(beta0,gamma0,sig20,nu0)

# system.time(
for (ii in 1:length(ns)) {
  system.time({
    n <- ns[ii]
    message("---- n = ", n, " ----")
    # covers.el <- rep(0,5)
    # lengths.el <- rep(0,5)
    # covers.hlm <- rep(0,5)
    # lengths.hlm <- rep(0,5)
    covers.el <- c(45, 42, 32, 27, 45)
    covers.hlm <- c(46, 42, 48, 24, 50)
    lengths.el <- c(14.84395, 14.83234, 16.27844, 22.78626, 28.60688)
    lenghts.hlm <- c(15.07314, 15.07314, 17.47164, 23.02124, 24.67563)
    
    reptimes <- 527
    for (jj in 78:reptimes) {
      message("** jj = ", jj, " **")
      # cat ("Press [enter] to continue\n")
      # line <- readline()
      if (jj %% 50 == 0) message("jj = ", jj)
      X <- cbind(rep(1,n),rnorm(n))
      Z <- matrix(rnorm(1*n),n,1)
      eps <- gen_eps(n, dist = "chisq", df = 5, tau = alpha)$eps
      yy <- c(X %*% beta0 + sqrt(sig20)*exp(Z %*% gamma0)*eps)
      # plot(yy~X[,2], cex=0.3)
      # random censoring
      cc <- rnorm(n,mean=2.5,sd=1)
      deltas <- yy<=cc
      y <- yy
      message("censored pct = ", sum(1-deltas)/n)
      if (sum(1-deltas)/n > 0.2 || sum(1-deltas)/n < 0.1) {
        reptimes <- reptimes+1
        next
      }
      y[as.logical(1-deltas)] <- cc[as.logical(1-deltas)]
      hlmout <- hlm(y, deltas, X, cbind(1,Z))
      message("hlm conv = ", hlmout$conv)
      if (!hlmout$conv) {
        reptimes <- reptimes+1
        next
      }
      
      beta.hlm <- hlmout$coef$beta
      gamma.hlm <- hlmout$coef$gamma[2]*0.5
      sig2.hlm <- exp(hlmout$coef$gamma[1])
      nu.hlm <- quantile((y-X %*% beta.hlm)*exp(-Z %*% gamma.hlm)/sqrt(sig2.hlm),alpha)
      
      ## el method
      betaInit <- beta.hlm
      gammaInit <- gamma.hlm
      sig2Init <- sig2.hlm
      nuInit <- nu.hlm
      # run adaptive MCMC here
      postout <- qrls_cens.post_adapt(y,X,Z,deltas,alpha,nsamples,nburn,
                                      betaInit,gammaInit,sig2Init,nuInit,
                                      mwgSd = mwgSds)
      theta_chain <- postout$theta_chain
      if (all(theta_chain == 0.0)) {
        message("MCMC not converged.")
        reptimes <- reptimes+1
        next
      }
      cat("paccept = ", postout$paccept, "\n")
      if (any(postout$paccept < 0.4) || any(postout$paccept > 0.5) ) {
        reptimes <- reptimes+1
        next
      }
      qts <- rep(NA,5)
      for (kk in 1:5) {
        qts <- quantile(theta_chain[kk,],c(0.025,0.975))
        # gamma range need to time 2 since the estimate for gamma/2
        # if (kk == 3) qts <- 2*qts
        if (theta0[kk] >= qts[1] && theta0[kk] <= qts[2]) covers.el[kk] <- covers.el[kk] + 1
        lengths.el[kk] <- lengths.el[kk] + (qts[2]-qts[1])
      }
      
      ## hlm method
      theta.boot <- qrls_cens.boot_R(y,X,Z,deltas,
                                     alpha,beta.hlm,gamma.hlm,sig2.hlm,nu.hlm,
                                     nboot=1000)
      CI.beta <- sapply(1:length(beta.hlm), function(ii) {
        qrls_cens.bootCI_R(theta.hat = beta.hlm[ii], theta.boot = theta.boot$beta.boot[,ii])
      })
      CI.gamma <- sapply(1:length(gamma.hlm), function(ii) {
        qrls_cens.bootCI_R(theta.hat = gamma.hlm[ii], theta.boot = theta.boot$gamma.boot[,ii])
      })
      CI.sig2 <- qrls_cens.bootCI_R(theta.hat = sig2.hlm, theta.boot = theta.boot$sig2.boot)
      CI.nu <- qrls_cens.bootCI_R(theta.hat = nu.hlm, theta.boot = theta.boot$nu.boot)
      for (mm in 1:2) {
        if (theta0[mm] >= CI.beta[1,mm] && theta0[mm] <= CI.beta[2,mm]) covers.hlm[mm] <- covers.hlm[mm] + 1
        lengths.hlm[mm] <- lengths.hlm[mm] + CI.beta[2]-CI.beta[1]
      }
      if (theta0[3] >= CI.gamma[1] && theta0[3] <= CI.gamma[2]) covers.hlm[3] <- covers.hlm[3] + 1
      lengths.hlm[3] <- lengths.hlm[3] + CI.gamma[2]-CI.gamma[1]
      if (theta0[4] >= CI.sig2[1] && theta0[4] <= CI.sig2[2]) covers.hlm[4] <- covers.hlm[4] + 1
      lengths.hlm[4] <- lengths.hlm[4] + CI.sig2[2]-CI.sig2[1]
      if (theta0[5] >= CI.nu[1] && theta0[5] <= CI.nu[2]) covers.hlm[5] <- covers.hlm[5] + 1
      lengths.hlm[5] <- lengths.hlm[5] + CI.nu[2]-CI.nu[1]
      
      message("reptimes = ", reptimes)
      cat("covers = ", covers.el, "\n")
      cat("lengths = ", lengths.el, "\n")
      cat("covers = ", covers.hlm, "\n")
      cat("lengths = ", lengths.hlm, "\n")
    }
    for (ll in 1:5) {
      coverProbs.el[ll,ii] <- covers[ll]/500
      lengthCI.el[ll,ii] <- lengths[ll]/500
    }
  })
}

# Restore output to console
sink() 
sink(type="message")

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


