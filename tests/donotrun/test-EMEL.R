#--------- check that EMEL.R is working properly -----------------
require(bayesEL)
source("~/bayesEL/tests/testthat/el-utils.R")
source("~/bayesEL/tests/donotrun/mle-check.R")

# 1-d problem --- still sometimes not stable
p <- 1
N <- 5
mu <- 1
z <- rnorm(N, mean = mu) # true lifetime variable
c <- rnorm(N, mean = 2*mu) # censoring variable 
delta <- z <= c
sum(!delta)/N # censored percentage
y <- z 
y[!delta] <- c[!delta] # observed lifetime
X <- matrix(rep(p,N),p,N) # p x n matrix
ws0 <- rep(1/N,N)
# order data here
ord <- order(y)
y <- y[ord]
X <- matrix(X[,ord],p,N)
delta <- delta[ord]
G <- mr.evalG_R(y, X, mu)
ws0 <- ws0[ord]

# optimization in C++
ws.cpp <- EMEL(G, delta, ws0)
ws.cpp
# optimization in R
wsout <- EMEL_R(G, delta, ws0)
ws.R <- wsout$ws
ws.R
# difference of the two implementation 
ws.cpp - ws.R


###### TODO: modify below ######

#---- 5 variables case: check the result of em_el ----
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
    wsps <- rep(NA,5)
    for (ii in 1:length(ws)) {
        wsps[ii] <- sum(ws[ii:5])
    }
    output <- list(ws=ws, 
                   objval=sum(delta*log(ws)+(1-delta)*log(wsps)))
    return(output)
}

# 2-dim problem 5 obvservations: 
p <- 2
N <- 5
beta <- c(2,5)
X <- rbind(rep(1,N),rnorm(N))
z <- beta %*% X + rnorm(N)
c <- rnorm(N, mean = 2*mean(z)) # censoring variable 
delta <- z <= c
sum(!delta)/N # censored percentage
y <- z 
y[!delta] <- c[!delta] # observed lifetime
ws0 <- rep(1/N,N)
# order data here
ord <- order(y)
y <- y[ord]
X <- matrix(X[,ord],p,N)
G <- mr.evalG_R(y, X, beta)
ws0 <- ws0[ord]

wsnew <- EMEL(G, delta, ws0)
wsnew_R <- EMEL_R(G, delta, ws0)$ws

# qs <- getqs_R(ws0,delta)

# m <- 2 # for mean regression always 2
# p <- 1
# n <- 5
# beta <- 1
# z <- rnorm(n, mean = mu)
# c <- rnorm(n, mean = 2*mu)
# delta <- z <= c
# sum(!delta)/n # censored percentage
# y <- z
# y[!delta] <- c[!delta]
# X <- matrix(rep(1,5),1,5) # p x n matrix
# G <- mr.evalG_R(y, X, mu)
# lambda0 <- rnorm(m)
# ws0 <- rep(1/5,5)
# # order data here
# ord <- order(y)
# y <- y[ord]
# X <- matrix(X[,ord],1,n)
# ws0 <- ws0[ord]

d <- seq(0,1,length.out = 100)
objval <- rep(NA, 100)
objvalmode <- -Inf
modeval <- NA
modews <- rep(NA,5)
for (ii in 1:100){
    wsobjval <- obj(d[ii],wsnew_R[5])
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
abline(v=wsnew[4], col='blue')
legend('topright',legend=c(expression('true mode'),
                           expression('emel result w[4]')),
       lty = c(1,1), col = c('red','blue'), cex = 0.4)


