# ---- experiments with adjusted G matrix ----
library(bayesEL)
library(optimCheck)
source("../testthat/el-utils.R")
source("../testthat/el-rfuns.R")
source("../testthat/el-model.R")
source("gen_eps.R")

# ---- support modification when there is no censoring ----
# mr with 1 intercept & 1 slope 
n <- 300
X0 <- matrix(rep(1,n),n,1)
X1 <- matrix(rnorm(n),n,1)
X <- cbind(X0,X1)
eps <- gen_eps(n)
beta_intercept <- 1
beta_slope <- 1.5
y <- 1 + c(X1 %*% beta_slope) + eps 
plot(X1,y,cex=0.3)
beta0 <- c(beta_intercept, beta_slope)

adjust <- TRUE
numpoints <- 100
beta1.seq <- seq(beta0[1]-3,beta0[1]+3,length.out = numpoints)
beta2.seq <- seq(beta0[2]-3,beta0[2]+3,length.out = numpoints)
beta.seq <- cbind(beta1.seq,beta2.seq)
logel.seq <- matrix(rep(NA,2*numpoints),2,numpoints)
for (ii in 1:numpoints) {
  if (ii %% 10 == 0) message("ii = ", ii)
  
  # beta0
  G <- mr.evalG(y,X,c(beta1.seq[ii],beta0[2]))
  if(adjust) {
    G <- adjG_R(G)
    omegas <- omega.hat_R(G,adjust=adjust)
    logel.seq[1,ii] <- logEL_R(omegas)
  }
  else {
    omegas <- omega.hat(G)
    logel.seq[1,ii] <- logEL(omegas)
  }
  
  # beta1
  G <- mr.evalG(y,X,c(beta0[1],beta2.seq[ii]))
  if(adjust) {
    G <- adjG_R(G)
    omegas <- omega.hat_R(G,adjust=adjust)
    logel.seq[2,ii] <- logEL_R(omegas,adjust)
  }
  else {
    omegas <- omega.hat(G)
    logel.seq[2,ii] <- logEL(omegas)
  }
}

par(mfrow=c(1,2))
if (!adjust) {
  # plot(logel.seq[1,])
  # plot(logel.seq[2,])
  logelmode1 <- plotEL(beta1.seq, logel.seq[1,], beta0[1], NA, 
                       log.scale = TRUE, expression(beta[0]), legend.loc='topright')
  logelmode2 <- plotEL(beta2.seq, logel.seq[2,], beta0[2], NA, 
                       log.scale = TRUE, expression(beta[1]), legend.loc='topright')
}
if (adjust) {
  # plot(logel.seq[1,])
  # plot(logel.seq[2,])
  logelmode1_adj <- plotEL(beta1.seq, logel.seq[1,], beta0[1], NA, 
                           log.scale = TRUE, expression(beta[0]), legend.loc='topright')
  logelmode2_adj <- plotEL(beta2.seq, logel.seq[2,], beta0[2], NA, 
                           log.scale = TRUE, expression(beta[1]), legend.loc='topright')
}


# ---- rand G ----
adjust <- TRUE
for (ii in 1:100) {
  if (ii %% 10 == 0) message("ii = ", ii)
  n <- 20
  p <- 10
  max_iter <- 500
  rel_tol <- 1e-7
  G <- matrix(rnorm(n*p), n, p)
  # omegahat.R <- omega.hat.NC_R(G, max_iter = max_iter, rel_tol = rel_tol, verbose = FALSE)
  
  if (adjust) G <- adjG_R(G)
  deltas <- rep(1,n)
  numcens <- sample(round(n/2),1)
  censinds <- sample(n,numcens)
  deltas[censinds] <- 0
  
  epsilons <- gen_eps(n)
  omegahat.R <- omega.hat.EM_R(G, deltas, epsilons, adjust = adjust,
                               max_iter = max_iter, rel_tol = rel_tol, verbose = FALSE)$omegas
  if (anyNA(omegahat.R)) break
}

# ---- use random G matrix to check if the logEL are close enough, number of steps to take 
dif <- c()
Gs <- list()
del <- list()
eps <- list()
nrep <- 100
count <- 0
while(TRUE) {
  count <- count + 1
  if (count == nrep+1) break
  message("count = ", count)
  n <- 100
  p <- 3
  max_iter <- 200
  rel_tol <- 1e-3
  G <- matrix(rnorm(n*p), n, p)
  deltas <- rep(1,n)
  numcens <- 15
  censinds <- sample(n,numcens)
  deltas[censinds] <- 0
  # if (sum(1-deltas)/n < 0.15 || sum(1-deltas)/n > 0.20) next
  epsilons <- gen_eps(n)
  omegas <- omega.hat.EM_R(G,deltas,epsilons)$omegas
  logel <- logEL(omegas,epsilons,deltas)
  G <- adjG_R(G)
  omegas <- omega.hat.EM_R(G, deltas, epsilons, adjust = TRUE,
                           max_iter = max_iter, rel_tol = rel_tol, verbose = FALSE)$omegas
  logel.adj <- logEL_R(omegas,epsilons,deltas,adjust=TRUE)
  
  if(!is.infinite(logel)) {
    Gs[[count]] <- G
    del[[count]] <- deltas
    eps[[count]] <- epsilons
    dif <- c(dif, logel-logel.adj)
  }
  else {
    count <- count-1
    next
  }
}

count
plot(dif,cex=.3)
which(abs(dif) > 1e-3)

idx <- 47
deltas <- del[[idx]]
epsilons <- eps[[idx]]
G <- Gs[[idx]][1:length(deltas),]

omegas <- omega.hat(G,deltas,epsilons)
logel <- logEL(omegas,epsilons,deltas)

G <- adjG_R(G)
omegas.adj <- omega.hat.EM_R(G, deltas, epsilons, adjust = TRUE,
                         max_iter = max_iter, rel_tol = rel_tol, verbose = FALSE)$omegas
logel.adj <- logEL_R(omegas.adj,epsilons,deltas,adjust=TRUE)

# if the artificial point got a large weight, then the likelihood is affected

# ---- mr model: 1-d case ----
n <- 100
mu0 <- rnorm(1,0,1)
# X <- matrix(rnorm(n,0,1),n,1) # each row of X is one observation
X <- matrix(rep(1,n), n, 1)

# dist is one of "norm","t","chisq","lnorm"
eps <- gen_eps(n, dist = "norm", df = NULL)

yy <- c(X * mu0) + eps
# random censoring
cc <- rnorm(n,mean=1,sd=1)
deltas <- yy<=cc
y <- yy
sum(1-deltas)/n
y[as.logical(1-deltas)] <- cc[as.logical(1-deltas)]

adjust <- TRUE
numpoints <- 100
mu.seq <- seq(-.5+mu0,.5+mu0,length.out = numpoints)
logel.seq <- rep(NA,numpoints)
for (ii in 1:numpoints) {
  if (ii %% 10 == 0) message("ii = ", ii)
  G <- mr.evalG(y,X,mu.seq[ii])
  epsilons <- y - c(X * mu.seq[ii])
  if (adjust) {
    G <- adjG_R(G)
    omegas <- omega.hat.EM_R(G, deltas, epsilons, adjust = adjust)$omegas
    logel.seq[ii] <- logEL_R(omegas,epsilons,deltas,adjust = adjust)
  }
  else {
    omegas <- omega.hat(G, deltas, epsilons)
    logel.seq[ii] <- logEL(omegas, epsilons, deltas)
  }
}
logelmode <- plotEL(mu.seq, logel.seq, mu0, mean(y), expression(mu))
# logel.seq.adj <- logel.seq

# ---- mr model: 2-d case ----
n <- 100
p <- 2
X1 <- matrix(rnorm(n),n,1)
X <- cbind(1,X1)

# dist is one of "norm","t","chisq","lnorm"
# eps <- gen_eps(n, dist = "lnorm", df = 5)
eps <- gen_eps(n)

beta_I <- 1
beta_S <- 1.5
# beta_I <- rnorm(1)
# beta_S <- rnorm(1)
beta0 <- c(beta_I, beta_S)
# plot(X1,y,cex=0.3)

# random censoring
# yy <- beta_I + c(X1 %*% beta_S) + eps 
# cc <- rnorm(n,mean=1.5*beta_I,sd=1)
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

# grid plot of conditionals: beta1|beta2 and beta2|beta1
adjust <- TRUE
numpoints <- 100
beta1.seq <- seq(beta0[1]-1,beta0[1]+1,length.out = numpoints)
beta2.seq <- seq(beta0[2]-1,beta0[2]+1,length.out = numpoints)
beta.seq <- cbind(beta1.seq,beta2.seq)
logel.seq <- matrix(rep(NA,2*numpoints),2,numpoints)
for (ii in 1:numpoints) {
  if (ii %% 10 == 0) message("ii = ", ii)
  G <- mr.evalG(y,X,c(beta1.seq[ii],beta0[2]))
  epsilons <- y - c(X %*% c(beta1.seq[ii],beta0[2]))
  if(adjust) {
    G <- adjG_R(G)
    # logel.seq[1,ii] <- logEL_EMAC_R(G,c(epsilons,-Inf),c(deltas,0))
    omegas <- omega.hat.EM_R(G,deltas,epsilons,adjust,dbg = TRUE,verbose = FALSE)$omegas
    logel.seq[1,ii] <- logEL_R(omegas,epsilons,deltas,adjust)
  }
  else {
    omegas <- omega.hat(G,deltas,epsilons)
    logel.seq[1,ii] <- logEL(omegas,epsilons,deltas)
  }
  
  G <- mr.evalG(y,X,c(beta0[1],beta2.seq[ii]))
  epsilons <- y - c(X %*% c(beta0[1],beta2.seq[ii]))
  if(adjust) {
    G <- adjG_R(G)
    # logel.seq[2,ii] <- logEL_EMAC_R(G,c(epsilons,-Inf),c(deltas,0))
    omegas <- omega.hat.EM_R(G,deltas,epsilons,adjust)$omegas
    logel.seq[2,ii] <- logEL_R(omegas,epsilons,deltas,adjust)
  }
  else {
    omegas <- omega.hat(G,deltas,epsilons)
    logel.seq[2,ii] <- logEL(omegas,epsilons,deltas)
  }
}

if (!adjust) {
  logelmode1 <- plotEL(beta1.seq, logel.seq[1,], beta0[1], NA, expression(beta[0]))
  logelmode2 <- plotEL(beta2.seq, logel.seq[2,], beta0[2], NA, expression(beta[1]))
}
if (adjust) {
  logelmode1_adj <- plotEL(beta1.seq, logel.seq[1,], beta0[1], NA, expression(beta[0]))
  logelmode2_adj <- plotEL(beta2.seq, logel.seq[2,], beta0[2], NA, expression(beta[1]))
}

any(is.infinite(logel.seq[2,]))

# --------------------------------------------------------------------------- 
# mr: is the conditional log EL curve smooth if there is no censoring for beta1?
n <- 100
p <- 2
X1 <- matrix(rnorm(n),n,1)
X <- cbind(1,X1)
beta_I <- 1
beta_S <- 1.5
beta0 <- c(beta_I, beta_S)
eps <- gen_eps(n)
y <- beta_I + c(X1 %*% beta_S) + eps 

adjust <- TRUE
numpoints <- 100
beta1.seq <- seq(beta0[1]-5,beta0[1]+5,length.out = numpoints)
beta2.seq <- seq(beta0[2]-5,beta0[2]+5,length.out = numpoints)
beta.seq <- cbind(beta1.seq,beta2.seq)
logel.seq <- matrix(rep(NA,2*numpoints),2,numpoints)

for (ii in 1:numpoints) {
  if (ii %% 10 == 0) message("ii = ", ii)
  G <- mr.evalG(y,X,c(beta1.seq[ii],beta0[2]))
  epsilons <- y - c(X %*% c(beta1.seq[ii],beta0[2]))
  if(adjust) {
    G <- adjG_R(G)
    omegas <- omega.hat_R(G, adjust=adjust)
    logel.seq[1,ii] <- logEL_R(omegas)
  }
  else {
    omegas <- omega.hat(G)
    logel.seq[1,ii] <- logEL(omegas)
  }
  
  G <- mr.evalG(y,X,c(beta0[1],beta2.seq[ii]))
  epsilons <- y - c(X %*% c(beta0[1],beta2.seq[ii]))
  if(adjust) {
    G <- adjG_R(G)
    omegas <- omega.hat_R(G, adjust=adjust)
    logel.seq[2,ii] <- logEL_R(omegas)
  }
  else {
    omegas <- omega.hat(G)
    logel.seq[2,ii] <- logEL(omegas)
  }
}

par(mfrow=c(1,2))
if (!adjust) {
  logelmode1 <- plotEL(beta1.seq, logel.seq[1,], beta0[1], NA, expression(beta[0]))
  logelmode2 <- plotEL(beta2.seq, logel.seq[2,], beta0[2], NA, expression(beta[1]))
}

# --------------------------------------------------------------------------- 
# mrls: is the conditional log EL curve smooth if there is no censoring for params other than beta0?
