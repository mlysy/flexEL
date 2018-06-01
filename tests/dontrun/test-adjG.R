# ---- experiments with adjusted G matrix ----
library(bayesEL)
source("../testthat/el-utils.R")
source("../testthat/el-rfuns.R")
source("../testthat/el-model.R")
source("gen_eps.R")

# ---- mr model 1-d case ----
for (ii in 1:1000) {
  n <- 20
  p <- 10
  max_iter <- 500
  rel_tol <- 1e-7
  adjust <- TRUE
  G <- matrix(rnorm(n*p), n, p)
  omegahat.R <- omega.hat.NC_R(G, max_iter = max_iter, rel_tol = rel_tol, verbose = FALSE)
  omegahat.R
  if (adjust) G <- adjG_R(G)
  deltas <- rep(1,n)
  numcens <- sample(round(n/2),1)
  censinds <- sample(n,numcens)
  deltas[censinds] <- 0
  sum(1-deltas)/n
  epsilons <- rnorm(n)

  omegahat.R <- omega.hat.EM_R(G, deltas, epsilons, adjust = adjust,
                               max_iter = max_iter, rel_tol = rel_tol, verbose = FALSE)
  if (anyNA(omegahat.R)) break
}

if (adjust) {
  n <- n+1
  epsilons <- c(epsilons,0)
  deltas <- c(deltas,0)
}
idx0 <- (abs(omegahat.R) < 1e-3 & !deltas)
ocheck <- optim_proj(xsol = rep(1,n-p-sum(idx0)),
                     xrng = 0.02,
                     npts = 201, 
                     fun = function(x) {omega.pcheck(x, omegahat.R, G, deltas, epsilons, idx0, rel_tol)},
                     plot = TRUE)


# generate data, compare the log likelihood with and without the adjusted point
n <- 200
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
# if (adjust) deltas <- c(deltas,0)
for (ii in 1:numpoints) {
  if (ii %% 10 == 0) message("ii = ", ii)
  G <- mr.evalG(y,X,mu.seq[ii])
  if (adjust) G <- adjG_R(G)
  epsilons <- y - c(X * mu.seq[ii])
  # if (adjust) epsilons <- c(epsilons,0)
  omegas <- omega.hat.EM_R(G, deltas, epsilons, adjust = adjust)
  logel.seq[ii] <- logEL_R(omegas,epsilons,deltas,adjust = adjust)
}
logelmode <- plotEL(mu.seq, logel.seq, mu0, mean(y), expression(mu))

adjust <- FALSE
for (ii in 1:numpoints) {
  if (ii %% 10 == 0) message("ii = ", ii)
  G <- mr.evalG(y,X,mu.seq[ii])
  if (adjust) G <- adjG_R(G)
  epsilons <- y - c(X * mu.seq[ii])
  # if (adjust) epsilons <- c(epsilons,0)
  omegas <- omega.hat.EM_R(G, deltas, epsilons, adjust = adjust)
  logel.seq[ii] <- logEL_R(omegas,epsilons,deltas,adjust = adjust)
}
logelmode <- plotEL(mu.seq, logel.seq, mu0, mean(y), expression(mu))

# ---- mr model: 2-d case ----

n <- 200
p <- 2
X1 <- matrix(rnorm(n),n,1)
X <- cbind(1,X1)
# eps <- rnorm(n) # N(0,1) error term

# dist is one of "norm","t","chisq","lnorm"
eps <- gen_eps(n, dist = "norm", df = NULL)

# beta_I <- 1
# beta_S <- 1.5
beta_I <- rnorm(1)
beta_S <- rnorm(1)
yy <- beta_I + c(X1 %*% beta_S) + eps 
beta0 <- c(beta_I, beta_S)
# plot(X1,y,cex=0.3)

# random censoring
cc <- rnorm(n,mean=2.5*beta_I,sd=1)
deltas <- yy<=cc
y <- yy
sum(1-deltas)/n
y[as.logical(1-deltas)] <- cc[as.logical(1-deltas)]

# grid plot of conditionals: beta1|beta2 and beta2|beta1
adjust <- FALSE
numpoints <- 100
beta1.seq <- seq(beta0[1]-.5,beta0[1]+.5,length.out = numpoints)
beta2.seq <- seq(beta0[2]-.5,beta0[2]+.5,length.out = numpoints)
beta.seq <- cbind(beta1.seq,beta2.seq)
logel.seq <- matrix(rep(NA,2*numpoints),2,numpoints)
for (ii in 1:numpoints) {
  if (ii %% 10 == 0) message("ii = ", ii)
  G <- mr.evalG(y,X,c(beta1.seq[ii],beta0[2]))
  if(adjust) G <- adjG_R(G)
  epsilons <- y - c(X %*% c(beta1.seq[ii],beta0[2]))
  omegas <- omega.hat.EM_R(G,deltas,epsilons,adjust)
  logel.seq[1,ii] <- logEL_R(omegas,epsilons,deltas,adjust)
  
  G <- mr.evalG(y,X,c(beta0[1],beta2.seq[ii]))
  if(adjust) G <- adjG_R(G)
  epsilons <- y - c(X %*% c(beta0[1],beta2.seq[ii]))
  omegas <- omega.hat.EM_R(G,deltas,epsilons,adjust)
  logel.seq[2,ii] <- logEL_R(omegas,epsilons,deltas,adjust)
}
logelmode1 <- plotEL(beta1.seq, logel.seq[1,], beta0[1], NA, expression(beta[0]))
logelmode2 <- plotEL(beta2.seq, logel.seq[2,], beta0[2], NA, expression(beta[1]))
