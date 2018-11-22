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

# 3-d plot with plotly
# to set labels: https://plot.ly/r/figure-labels/
# type could be "surface" or "contour"
plot3d <- function(type="surface", seq.x, seq.y, seq.z, 
                   xlab="intercept", ylab="slope", zlab="logEL",
                   title=NULL,
                   m1, logel.m1, lab.m1="true val",
                   m2, logel.m2, lab.m2="est val") {
  p <- plot_ly(x=seq.x, y=seq.y, z=seq.z,
               type=type) %>% 
    layout(
      title = title,
      scene = list(
        xaxis = list(title = xlab),
        yaxis = list(title = ylab),
        zaxis = list(title = zlab))
      )
  if (!missing(m1) && !missing(m2)) {
    p %>% add_trace(x=m1[1], y=m1[2], z=logel.m1,
                    name = "true val",
                    type = "scatter3d",
                    mode = "markers",
                    marker = list(color = "red")) %>% 
      add_trace(x=m2[1], y=m2[2], z=logel.m2,
                name = "est val",
                type = "scatter3d",
                mode = "markers",
                marker = list(color = "blue"))
  }
  else p
}

# EL curve plot
plotEL <- function(mu.seq, logel.seq, trueval, obs = NA, mu.name = 'param', legend.loc = 'topright') {
    plot(mu.seq, exp(logel.seq-max(logel.seq)),
         cex=0.2, xlab = mu.name, ylab = 'log EL', type = 'l')
    abline(v = trueval, col = 'red', lty=2) # true param
    abline(v = mu.seq[which.max(logel.seq)], lty=2) # mode of EL
    if (!is.na(obs)) {
        abline(v = obs, col='blue', lty=2) # observed mean / quantile
        legend('topright',legend=c('true param', 'logEL mode', 'observed'),
               lty = c(2,2,2), col = c('red','black','blue'), cex = 0.6)
    }
    else{
        legend(legend.loc,legend=c('true param', 'logEL mode'),
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
    
    logel <- rep(NA,n)
    logel[deltas==1] <- log(xNG[deltas==1])
    logel[deltas==0] <- log(psos[deltas==0])
    return(sum(logel))
    
    # return(sum(deltas*log(xNG)+(1-deltas)*log(psos)))
  }
}

# partial check: only on the omegas that are significantly nonzero
omega.pcheck <- function(x, omegas, G, deltas, epsilons, idx0, rel_tol) {
  p <- ncol(G)
  n <- nrow(G)
  G <- G[!idx0,]
  NG <- Null(G)
  xNG <- c(NG %*% x) - rowSums(NG) + omegas[!idx0] # x == 1s, xNG == omegas 
  if (any(xNG < -1e-10)) return(-Inf)
  xNG <- abs(xNG) # try replace the extreme small value to positive
  xNG <- xNG / sum(xNG) # normalize it
  if (missing(deltas) && missing(epsilons)) { # actually this should not be needed
    return(sum(log(xNG)))
  }
  else {
    xNGfull <- rep(0,n)
    xNGfull[!idx0] <- xNG
    epsOrd <- order(epsilons) # ascending order of epsilons
    # print(epsOrd)
    psos <- rep(0,n)
    for (ii in 1:n) {
      psos[ii] <- evalPsos_R(ii, epsOrd, xNGfull) 
    }
    logel <- rep(NA,n)
    logel[deltas==1] <- log(xNGfull[deltas==1])
    logel[deltas==0] <- log(psos[deltas==0])
    return(sum(logel))
  }
}