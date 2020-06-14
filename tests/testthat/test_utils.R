
# max of min of abs and rel error
max_xdiff <- function(x) {
  xdiff <- abs(diff(x))
  max(pmin(xdiff[,1], xdiff[,2]))
}

# library(MASS) # for use of Null
# wrapper of omega_hat_R for optimCheck
# if omegas is optimal, then x == rep(1,n-p) should be optimal
omega_check <- function(x, omegas, G, deltas, epsilons) {
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
omega_pcheck <- function(x, omegas, G, deltas, epsilons, idx0, rel_tol) {
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

# library(MASS) # for use of Null
# wrapper of omega_hat_R for optimCheck
# if omegas is optimal, then x == rep(1,n-p) should be optimal
omega_smooth_check <- function(x, omegas, G, deltas, epsilons, s=10) {
  NG <- Null(G)
  xNG <- c(NG %*% x) - rowSums(NG) + omegas # x == 1s, xNG == omegas 
  # might have to use a small negative value to aviod rounding errors
  if (any(xNG < -1e-7)) return(-Inf)
  xNG <- abs(xNG) # try replace the extreme small value to positive
  xNG <- xNG / sum(xNG) # normalize it
  return(logEL_smooth_R(xNG,epsilons,deltas,s))
}

omega_smooth_pcheck <- function(x, omegas, G, deltas, epsilons,idx0, s=10) {
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
    # xNGfull[idx0] <- omegas[idx0]
    # return(logEL_smooth_R(xNGfull,epsilons,deltas,s))
    # epsOrd <- order(epsilons) # ascending order of epsilons
    psos <- rep(0,n)
    for (ii in 1:n) {
      psos[ii] <- evalPsos_smooth_R(ii,xNGfull,epsilons,s)
    }
    logel <- rep(NA,n)
    logel[deltas==1] <- log(xNGfull[deltas==1])
    logel[deltas==0] <- log(psos[deltas==0])
    return(sum(logel))
  }
}