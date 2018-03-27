## check that omega.hat is working properly in low dimension cases 
require(bayesEL)
require(optimCheck)
source("~/bayesEL/tests/testthat/el-utils.R")
# source("~/bayesEL/tests/testthat/mle-check.R")

#---- 1-d problem -- still sometimes not stable ----
p <- 1
N <- 5
mu0 <- 1
z <- rnorm(N, mean = mu0) # true lifetime variable
c <- rnorm(N, mean = 2*mu0) # censoring variable 
deltas <- z <= c
sum(!deltas)/N # censored percentage
y <- z 
y[!deltas] <- c[!deltas] # observed lifetime
X <- matrix(rep(p,N),N,p) # p x n matrix
mu <- mu0+rnorm(1)
epsilons <- y - mu
G <- mr.evalG_R(y, X, mu)

# optimization in R
omegas.R <- omega.hat_R(G, deltas, epsilons)
omegas.R

# optimization in C++
omegas.cpp <- omega.hat(G, deltas, epsilons)
omegas.cpp

# difference of the two implementation 
omegas.cpp - omegas.R

# So .. if they both converged then they are the same but often only one converged
# but, how to check optimality .. 

# ok choose which to check here
omegas <- omegas.cpp
# quick check by using Null space vectors
NG <- Null(G) # any col c of NG satisfies c %*% G = 0
# we know that omegas > 0, just take any c of NG, add omegas until > 0
# then normalize it
ind <- sample(N-p,1) # dim of NG is n x (n-p)
mult <- sample(20:100,1) # randomly scale up the omegas by mult
omegas_twk <- NG[,ind] + mult*omegas
# make sure omegas_twk is nonnegative: 
while (sum(omegas_twk < 0) > 0) omegas_twk <- omegas_twk + 50*omegas
omegas_twk <- omegas_twk / sum(omegas_twk) # normalize it
logel <- logEL_R(omegas, G, deltas, epsilons)
logel_twk <- logEL_R(omegas_twk, G, deltas, epsilons)
logel - logel_twk
# expect_gt(logel.cpp, logel.cpp_twk)

#---- 2-d priblem -- 5 variables case ----
# Note: epsilons should be in the enviornment 
obj  <- function(d,e) {
    G <- t(G)
    c2 <- G[1,1]-G[1,2]
    c3 <- G[1,1]-G[1,3]
    c4 <- G[1,4]-G[1,1]
    c5 <- G[1,5]-G[1,1]
    d2 <- G[1,2]*G[2,1]-G[1,1]*G[2,2]
    d3 <- G[1,3]*G[2,1]-G[1,1]*G[2,3]
    d4 <- G[1,1]*G[2,4]-G[1,4]*G[2,1]
    d5 <- G[1,1]*G[2,5]-G[1,5]*G[2,1]
    b <- 1/(c2*d3-c3*d2)*(G[1,1]*d3 + (c4*d3-c3*d4)*d + (c5*d3-c3*d5)*e)
    c <- 1/(c3*d2-c2*d3)*(G[1,1]*d2 + (c4*d2-c2*d4)*d + (c5*d2-c2*d5)*e)
    a <- 1-(b+c+d+e)
    ws <- c(a,b,c,d,e)
    # DEBUG
    # if(sum(ws < 0) > 0) {
    #   message("negative weights")
    # }
    # else {
    #   message("ok")
    # }
    # sum omegas in the order of epsilons
    weights <- evalWeights_R(deltas, ws, epsilons)
    objval <- sum(weights * log(ws))
    output <- list(ws=ws, 
                   objval=objval)
    return(output)
}

# simulate data
p <- 2
N <- 5
beta0 <- c(5,2)
X <- cbind(rep(1,N),rnorm(N)) # nObs x nEqs
z <- c(X %*% beta0 + rnorm(N))
c <- rnorm(N, mean = 1.5*mean(z)) # censoring variable
deltas <- z <= c
deltas
# sum(!deltas)/N # censored percentage
y <- z 
y[!deltas] <- c[!deltas] # observed lifetime
beta <- beta0 + rnorm(p)
epsilons <- c(y - X %*% beta)
G <- mr.evalG_R(y, X, beta)

omegas.R <- omega.hat_R(G, deltas, epsilons)
omegas.R

omegas.cpp <- omega.hat(G, deltas, epsilons)
omegas.cpp

# difference of the two 
omegas.cpp - omegas.R 

d <- seq(0,1,length.out = 100)
objval <- rep(NA, 100)
objvalmode <- -Inf
modeval <- NA
modews <- rep(NA,5)
## choose which to plot here 
omegas_plot <- omegas.cpp
for (ii in 1:100){
    wsobjval <- obj(d[ii], 0.05) # omegas_plot[5]
    objval[ii] <- wsobjval$objval
    # if(!is.na(objval[ii])) print(objval)
    if (!is.na(objval[ii]) && (objval[ii] > objvalmode)) {
        objvalmode <- objval[ii]
        modeval <- d[ii]
        modews <- wsobjval$ws
    }
}
objval
plot(d,objval,type = 'l')
abline(v=modeval, col='red')
abline(v=omegas_plot[4], col='blue')
legend('topright',
       legend=c(expression('true mode'),
                expression('emel result w[4]')),
       lty = c(1,1), col = c('red','blue'), cex = 0.4)


