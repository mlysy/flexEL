# see if this smoothes out the censored log likelihood 

ind <- function(x) {
  return(as.numeric(x <= 0))
}

# smoothed indicator function 1(x <= 0)
ind.smooth_R <- function(x, s=10) {
  return(1/(1+exp(s*x)))
}

# 1st derivative of smoothed indicator function 1(x <= 0)
ind1.smooth_R <- function(x, s=10) {
  return(-s*exp(s*x)/(1+exp(s*x))^2)
}

# sfun <- function(x,s) {
#   ind1.smooth_R(x=s,s=x)
# }

# smoothed partial sum of omegas (for loop)
# evalPsos.smooth_R <- function(ii, omegas, epsilons, s=10) {
#   # ei <- epsilons[ii]
#   n <- length(omegas)
#   psos <- 0
#   for (jj in 1:n) {
#     # ej <- epsilons[jj]
#     psos <- psos + ind.smooth_R(epsilons[ii]-epsilons[jj],s)*omegas[jj]
#   }
#   return(psos)
# }

# smoothed partial sum of omegas (vectorized version)
evalPsos.smooth_R <- function(ii, omegas, epsilons, s=10) {
  psos <- sum(ind.smooth_R(epsilons[ii]-epsilons,s)*omegas)
  return(psos)
}

# smoothed censored logEL
logEL.smooth_R <- function(omegas,epsilons,deltas,s=10) {
  epsOrd <- order(epsilons) # ascending order of epsilons
  n <- length(omegas)
  psos <- rep(0,n)
  for (ii in 1:n) {
    psos[ii] <- evalPsos.smooth_R(ii,omegas,epsilons,s) 
  }
  # numerical stability: watch out for extremely small negative values
  omegas[abs(omegas) < 1e-10/length(omegas)] <- 1e-10
  return(sum(deltas*log(omegas)+(1-deltas)*log(psos)))
}

# smoothed evalWeight function with sigmoid function (double for-loop)
# evalWeights.smooth_R <- function(deltas, omegas, epsilons, s=10) {
#   nObs <- length(omegas)
#   epsOrd <- order(epsilons)
#   psots <- rep(0,nObs)
#   for (ii in 1:nObs) {
#     for (jj in 1:nObs) {
#       psots[ii] <- psots[ii] + (1-deltas[jj])*ind.smooth_R(epsilons[jj]-epsilons[ii],s)*
#                                       omegas[ii]/evalPsos.smooth_R(jj,omegas,epsilons,s)
#     }
#   }
#   # the weights have the order as the original sample
#   weights <- deltas + psots
#   return(weights)
# }

# smoothed evalWeight function with sigmoid function (more vectorized version)
evalWeights.smooth_R <- function(deltas, omegas, epsilons, s=10) {
  nObs <- length(omegas)
  epsOrd <- order(epsilons)
  psots <- rep(0,nObs)
  psoss <- rep(NA,nObs)
  for (jj in 1:nObs) {
    psoss[jj] <- evalPsos.smooth_R(jj,omegas,epsilons,s)
  }
  for (ii in 1:nObs) {
    psots[ii] <- sum((1-deltas)*ind.smooth_R(epsilons-epsilons[ii],s)*omegas[ii]/psoss)
  }
  weights <- deltas + psots
  return(weights)
}

# using smoothed objective function in EM
omega.hat.EM.smooth_R <- function(G, deltas, epsilons, s=10, adjust = FALSE, 
                                  max_iter = 200, rel_tol = 1e-7, abs_tol = 1e-3, 
                                  verbose=FALSE, dbg = FALSE) {
  n <- nrow(G)
  m <- ncol(G)
  err <- Inf
  nIter <- 0
  # initialize omegas with uncensored solution 
  omegas <- omega.hat.NC_R(G, adjust, max_iter, rel_tol, verbose=FALSE)
  if (any(is.nan(omegas))) {
    message("Initial omegas are nans.")
    # return(rep(NaN,length(deltas)))
    return(list(conv=FALSE, omegas=rep(NaN,length(deltas))))
  }
  if (adjust) {
    epsilons <- c(epsilons,-Inf)
    deltas <- c(deltas,0)
  }
  logelOld <- logEL.smooth_R(omegas,epsilons,deltas,s)
  
  # for debug
  if (dbg) {
    omegamMat <- matrix(NA,nrow=n,ncol=max_iter+1)
    omegamMat[,1] <- omegas
    logels <- c(logelOld)
    artome <- c()
  }
  
  for (ii in 1:max_iter) {
    nIter <- ii
    # E step: calculating weights
    weights <- evalWeights.smooth_R(deltas, omegas, epsilons, s)
    # M step:
    # lambdaOut <- lambdaNRC_R(G, weights, max_iter, rel_tol, verbose=FALSE)
    lambdaOut <- lambdaNRC_R(G, weights, max_iter, rel_tol, verbose=FALSE)
    # TODO: what if not converged ?? use a random weights and continue ?
    if (!lambdaOut$convergence) {
      # message("lambdaNRC did not converge in EM")
      if (dbg) {
        return(list(conv=FALSE, omegas=rep(NaN,n), logels=logels, artome=artome))
      }
      else {
        return(list(conv=FALSE, omegas=rep(NaN,n)))
      }
    }
    lambdaNew <- lambdaOut$lambda
    qlg <- c(sum(weights) + lambdaNew %*% t(G))
    omegas <- weights/qlg
    # omegas <- omegas/sum(omegas)
    if (any(omegas < -rel_tol)) message("omega.hat.EM_R: negative omegas.")
    omegas <- abs(omegas)
    omegas <- omegas/sum(omegas)
    logel <- logEL.smooth_R(omegas,epsilons,deltas,s)
    
    # for debug
    if (dbg) {
      omegamMat[,ii+1] <- omegas
      logels <- c(logels,logel)
      artome <- c(artome,weights[n])
    }
    err <- abs(logel-logelOld)
    # if (verbose && nIter %% 20 == 0) {
    if (verbose) {
      message("nIter = ", nIter)
      message("abs err = ", abs(logel-logelOld))
    }
    if (err < abs_tol) break
    logelOld <- logel
  }
  notconv <- (nIter == max_iter && err > abs_tol) # TRUE if not converged 
  if (notconv) {
    # omegas = rep(NaN,n)
    if (dbg) {
      return(list(conv=FALSE,omegas = rep(NaN,n),logels=logels,artome=artome))
    }
    else {
      return(list(conv=FALSE,omegas = rep(NaN,n)))
    }
  }
  
  # for debug
  if (dbg) {
    return(list(conv=TRUE,
                nIter=nIter,
                omegas=omegas,
                logels=logels,
                omegamMat=omegamMat,
                artome=artome))
  }
  # else return(omegas)
  else {
    return(list(conv=TRUE, omegas=omegas, lambda=lambdaNew, weights=weights))
  }
}

# curve(ind.smooth_R(x,s=10),from=-5,to=5,col='red')
# curve(ind.smooth_R(x,s=100),from=-5,to=5,col='blue',add=TRUE)
# curve(ind(x),from=-5,to=5,col='black',add=TRUE)
# legend('topright', legend=c("indicator","smoothed s=10", "smoothed s = 100"), lty=c(1,1,1), col=c("black","red","blue"),cex=.8)

library(MASS) # for use of Null
# wrapper of omega.hat_R for optimCheck
# if omegas is optimal, then x == rep(1,n-p) should be optimal
omega.smooth.check <- function(x, omegas, G, deltas, epsilons, s=10) {
  NG <- Null(G)
  xNG <- c(NG %*% x) - rowSums(NG) + omegas # x == 1s, xNG == omegas 
  # might have to use a small negative value to aviod rounding errors
  if (any(xNG < -1e-7)) return(-Inf)
  xNG <- abs(xNG) # try replace the extreme small value to positive
  xNG <- xNG / sum(xNG) # normalize it
  return(logEL.smooth_R(xNG,epsilons,deltas,s))
}

omega.smooth.pcheck <- function(x, omegas, G, deltas, epsilons,idx0, s=10) {
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
    # return(logEL.smooth_R(xNGfull,epsilons,deltas,s))
    # epsOrd <- order(epsilons) # ascending order of epsilons
    psos <- rep(0,n)
    for (ii in 1:n) {
      psos[ii] <- evalPsos.smooth_R(ii,xNGfull,epsilons,s)
    }
    logel <- rep(NA,n)
    logel[deltas==1] <- log(xNGfull[deltas==1])
    logel[deltas==0] <- log(psos[deltas==0])
    return(sum(logel))
  }
}

# # sandwich estimator for convariance matrix
# library(numDeriv)
# # f_{ii} for mr_cens
# # ss is s here, just s seems to be a param in grad so changed the name
# mr_cens.logELii.smooth_R <- function(theta, y, X, deltas, ii, ss=10) {
#   nObs <- nrow(X)
#   nEqs <- ncol(X)
#   G <- mr.evalG(y,X,theta)
#   epsilons <- evalEpsilons_R(y,X,theta)
#   oout <- omega.hat.EM.smooth_R(G,deltas,epsilons,ss)
#   qs <- oout$weights
#   lambda <- oout$lambda
#   if (deltas[ii]) {
#     return(log(qs[ii]/(n+G[ii,] %*% lambda)))
#   }
#   else {
#     denom <- c(n + G %*% lambda)
#     return(sum(ind.smooth_R(epsilons[ii]-epsilons,ss)*qs/(denom)))
#   }
# }
# 
# mr_cens.sandCov.smooth_R <- function(y, X, deltas, beta.hat, s=10) {
#   nObs <- nrow(X)
#   nEqs <- ncol(X)
#   A <- hessian(function(b) {-mr_cens.neglogEL.smooth_R(y, X, deltas, b, s)}, x=beta.hat)
#   B <- matrix(0,nrow=nEqs,ncol=nEqs)
#   for (ii in 1:nObs) {
#     gii <- grad(mr_cens.logELii.smooth_R, x=beta.hat, y=y, X=X, deltas=deltas, ii=ii, ss=s)
#     B <- B + tcrossprod(gii,gii)
#   }
#   Ainv <- solve(A)
#   return(Ainv %*% B %*% Ainv)
# }

mr_cens.sandCov.smooth_R <- function(y, X, deltas, beta.hat, s=10) {
  nObs <- nrow(X)
  nEqs <- ncol(X)
  A <- hessian(function(b) {-mr_cens.neglogEL.smooth(y, X, deltas, b, s)}, x=beta.hat)
  # TODO: resample data for A as well or only for B?
  inds <- sample(nObs,replace = TRUE)
  X <- X[inds,]
  y <- y[inds]
  deltas <- deltas[inds]
  B <- grad(mr_cens.neglogEL.smooth, x=beta.hat, y=y, X=X, delta=deltas, sp=s)
  B <- tcrossprod(B,B)
  Ainv <- solve(A)
  return(Ainv %*% B %*% Ainv)
}


