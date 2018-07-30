# see if this smoothes out the censored log likelihood 

# smoothed indicator function 1(x <= 0)
ind.smooth_R <- function(x, s=10) {
  return(1/(1+exp(s*x)))
}

# smoothed partial sum of omegas
evalPsos.smooth_R <- function(ii, omegas, epsilons, s=10) {
  # ei <- epsilons[ii]
  n <- length(omegas)
  psos <- 0
  for (jj in 1:n) {
    # ej <- epsilons[jj]
    psos <- psos + ind.smooth_R(epsilons[ii]-epsilons[jj],s)*omegas[jj]
  }
  return(psos)
}

# smoothed censored logEL
logEL.smooth_R <- function(omegas,epsilons,deltas,s=10) {
  epsOrd <- order(epsilons) # ascending order of epsilons
  n <- length(omegas)
  psos <- rep(0,n)
  for (ii in 1:n) {
    psos[ii] <- evalPsos.smooth_R(ii, omegas, epsilons) 
  }
  # numerical stability: watch out for extremely small negative values
  omegas[abs(omegas) < 1e-10/length(omegas)] <- 1e-10
  return(sum(deltas*log(omegas)+(1-deltas)*log(psos)))
}

# smoothed evalWeight function with sigmoid function
evalWeights.smooth_R <- function(deltas, omegas, epsilons, s=10) {
  nObs <- length(omegas)
  epsOrd <- order(epsilons)
  psots <- rep(0,nObs)
  for (ii in 1:nObs) {
    for (jj in 1:nObs) {
      psots[ii] <- psots[ii] + (1-deltas[jj])*ind.smooth_R(epsilons[jj]-epsilons[ii],s)*
                                      omegas[ii]/evalPsos.smooth_R(jj,omegas,epsilons,s)
    }
  }
  # the weights have the order as the original sample 
  weights <- deltas + psots
  return(weights)
}

# using smoothed objective function in EM
omega.hat.EM.smooth_R <- function(G, deltas, epsilons, s=10, adjust = FALSE, 
                                  max_iter = 200, rel_tol = 1e-3, verbose=FALSE,
                                  dbg = FALSE) {
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
  # omegasOld <- omegas
  # logelOld <- logEL_R(omegas,epsilons,deltas)
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
    # err <- MaxRelErr(lambdaNew,lambdaOld)
    # err <- MaxRelErr(omegas,omegasOld)
    # logel <- logEL_R(omegas,epsilons,deltas)
    logel <- logEL.smooth_R(omegas,epsilons,deltas,s)
    
    # for debug
    if (dbg) {
      omegamMat[,ii+1] <- omegas
      logels <- c(logels,logel)
      artome <- c(artome,weights[n])
    }
    
    # err <- MaxRelErr(logel,logelOld)
    err <- abs(logel-logelOld)
    # if (verbose && nIter %% 20 == 0) {
    if (verbose) {
      message("nIter = ", nIter)
      # message("err = ", err)
      message("abs err = ", abs(logel-logelOld))
    }
    if (err < rel_tol) break
    logelOld <- logel
  }
  notconv <- (nIter == max_iter && err > rel_tol) # TRUE if not converged 
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