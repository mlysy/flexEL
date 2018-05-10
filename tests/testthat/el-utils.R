# source("el-rfuns.R")
# source("el-model.R")

#---- general helper functions ----

# generate matrix of iid N(0,1)
rMNorm <- function(n, p) {
  if(missing(p)) p <- n
  matrix(rnorm(n*p), n, p)
}

# max of min of abs and rel error
max.xdiff <- function(x) {
    xdiff <- abs(diff(x))
    max(pmin(xdiff[,1], xdiff[,2]))
}

plotEL <- function(mu.seq, logel.seq, trueval, meanobs = NA, mu.name = "param") {
    plot(mu.seq, exp(logel.seq-max(logel.seq)),
         cex=0.2, xlab = mu.name, ylab = 'log EL', type = 'l')
    abline(v = trueval, col = 'red', lty=2) # true param
    abline(v = mu.seq[which.max(logel.seq)], lty=2) # mode of EL
    if (!is.na(meanobs)) {
        abline(v = meanobs, col='blue', lty=2) # mean of observation
        legend('topleft',legend=c('true param', 'logEL mode', 'observed mean'),
               lty = c(2,2,2), col = c('red','black','blue'), cex = 0.6)
    }
    else{
        legend('topleft',legend=c('true param', 'logEL mode'),
               lty = c(2,2), col = c('red','black'), cex = 0.6)
    }
    return(mu.seq[which.max(logel.seq)]) # return the mode
}

norm_pdf <- function(ly, x) {
    py <- exp(ly - max(ly))
    py/sum(py)/(x[2]-x[1])
}

library(MASS) # for use of Null
# wrapper of omega.hat_R for optimCheck
# if omegas is optimal, then x == rep(1,n-p) should be optimal
omega.check <- function(x, omegas, G, deltas, epsilons) {
  NG <- Null(G)
  xNG <- c(NG %*% x) - rowSums(NG) + omegas # x == 1s, xNG == omegas 
  # might have to use a small negative value to aviod rounding errors
  if (any(xNG < -1e-10)) return(-Inf)
  xNG <- abs(xNG) # try replace the extreme small value to positive
  xNG <- xNG / sum(xNG) # normalize it
  if (missing(deltas) && missing(epsilons)) {
    return(sum(log(xNG)))
  }
  else {
    n <- length(omegas)
    epsOrd <- order(epsilons) # ascending order of epsilons
    # print(epsOrd)
    psos <- rep(0,n)
    for (ii in 1:n) {
      psos[ii] <- evalPsos_R(ii, epsOrd, xNG) 
    }
    # print(psos)
    return(sum(deltas*log(xNG)+(1-deltas)*log(psos)))
  }
}

