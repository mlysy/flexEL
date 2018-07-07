
# Aitken Acceleration
logEL_EMAC_R <- function(G, epsilons, deltas, 
                         max_iter = 100, rel_tol = 1e-3, verbose=FALSE, dbg=FALSE) {
  n <- nrow(G)
  m <- ncol(G)
  err <- Inf
  nIter <- 0
  # initialize omegas with uncensored solution 
  omegas <- omega.hat.NC_R(G, adjust, max_iter, rel_tol, verbose=FALSE)
  if (any(is.nan(omegas))) {
    message("Initial omegas are nans.")
    return(-Inf)
  }
  
  logelOld <- logEL_R(omegas,epsilons,deltas)
  logeltuple <- rep(0.0,3)
  
  for (ii in 1:max_iter) {
    nIter <- ii
    # E step: calculating weights
    weights <- evalWeights_R(deltas, omegas, epsilons)
    # M step:
    lambdaOut <- lambdaNRC_R(G, weights, max_iter, rel_tol, verbose=FALSE)
    # TODO: what if not converged ?
    if (!lambdaOut$convergence) {
      message("lambdaNRC did not converge in EM")
      return(-Inf)
    }
    lambdaNew <- lambdaOut$lambda
    qlg <- c(sum(weights) + lambdaNew %*% t(G))
    omegas <- weights/qlg
    omegas <- omegas/sum(omegas)
    if (any(omegas < -rel_tol)) message("omega.hat.EM_R: negative omegas.")
    omegas <- abs(omegas)
    
    logel <- logEL_R(omegas,epsilons,deltas)
    logeltuple[1:2] <- logeltuple[2:3]
    logeltuple[3] <- logel
    
    # TODO: use absolute error instead of relative error
    if (ii <= 3) {
      # err <- MaxRelErr(logel,logelOld)
      err <- abs(logel-logelOld)
      if (err < rel_tol) {
        epsOrd <- order(epsilons) # ascending order of epsilons
        n <- length(omegas)
        psos <- rep(0,n)
        for (ii in 1:n) {
          psos[ii] <- evalPsos_R(ii, epsOrd, omegas) 
        }
        # numerical stability: watch out for extremely small negative values
        omegas[abs(omegas) < 1e-10/length(omegas)] <- 1e-10
        message("Not accelerated.")
        if (verbose) {
          message("nIter = ", ii)
          message("err = ", err)
        }
        return(sum(deltas*log(omegas)+(1-deltas)*log(psos)))
      }
      else {
        logelOld <- logel
        if (ii == 3) {
          logelOld.ac <- logelOld
          if (dbg) logels <- c(logelOld.ac)
        }
      }
    }
    else {
      dif <- diff(logeltuple)
      c <- dif[2]/dif[1]
      # message("c = ", c)
      logel.ac  <- logeltuple[1]+1/(1-c)*dif[1]
      if (dbg) logels <- c(logels,logel.ac)
      # err <- MaxRelErr(logel.ac,logelOld.ac)
      err <- abs(logel.ac-logelOld.ac)
      if (verbose) {
        message("nIter = ", ii)
        # message("err = ", err)
        message("abs err = ", abs(logelOld.ac-logel.ac))
      }
      if (err < rel_tol) {
        if (dbg) return(list(logel=logel.ac, logels=logels))
        else return(logel.ac)
      }
      else {
        logelOld.ac <- logel.ac
      }
    }
    # if (verbose && nIter %% 20 == 0) {
  }
  if (ii == max_iter) return(-Inf)
}


# fitted acceleration
# logEL_EMAC_R <- function(G, epsilons, deltas, 
#                          max_iter = 100, rel_tol = 1e-5, fit_iter = 7, verbose=FALSE) {
#   n <- nrow(G)
#   m <- ncol(G)
#   err <- Inf
#   nIter <- 0
#   # initialize omegas with uncensored solution 
#   omegas <- omega.hat.NC_R(G, adjust, max_iter, rel_tol, verbose=FALSE)
#   if (any(is.nan(omegas))) {
#     message("Initial omegas are nans.")
#     return(-Inf)
#   }
#   
#   logelOld <- logEL_R(omegas,epsilons,deltas)
#   logels <- c(logelOld)
#   
#   for (ii in 1:max_iter) {
#     nIter <- ii
#     # E step: calculating weights
#     weights <- evalWeights_R(deltas, omegas, epsilons)
#     # M step:
#     lambdaOut <- lambdaNRC_R(G, weights, max_iter, rel_tol, verbose)
#     # TODO: what if not converged ?
#     if (!lambdaOut$convergence) {
#       message("lambdaNRC did not converge in EM")
#       return(-Inf)
#     }
#     lambdaNew <- lambdaOut$lambda
#     qlg <- c(sum(weights) + lambdaNew %*% t(G))
#     omegas <- weights/qlg
#     omegas <- omegas/sum(omegas)
#     if (any(omegas < -rel_tol)) message("omega.hat.EM_R: negative omegas.")
#     omegas <- abs(omegas)
#     logel <- logEL_R(omegas,epsilons,deltas)
#     logels <- c(logels,logel)
#     err <- MaxRelErr(logel,logelOld)
#     
#     if (err < rel_tol) {
#       epsOrd <- order(epsilons) # ascending order of epsilons
#       n <- length(omegas)
#       psos <- rep(0,n)
#       for (ii in 1:n) {
#         psos[ii] <- evalPsos_R(ii, epsOrd, omegas) 
#       }
#       # numerical stability: watch out for extremely small negative values
#       omegas[abs(omegas) < 1e-10/length(omegas)] <- 1e-10
#       message("Not accelerated.")
#       return(sum(deltas*log(omegas)+(1-deltas)*log(psos)))
#     }
#     else if (ii >= fit_iter) {
#       # browser()
#       idx <- 1/(1:(ii+1))
#       fit <- lm(logels ~ idx + I(idx^2))
#       err <- MaxRelErr(logels[(ii-2):(ii+1)],predict(fit)[(ii-2):(ii+1)])
#       if (err < 10*rel_tol) {
#         message("fit_iter = ", ii)
#         return(unname(coef(fit)[1]))
#       }
#     }
#     logelOld <- logel
#     # if (verbose && nIter %% 20 == 0) {
#     if (verbose) {
#       message("nIter = ", nIter)
#       message("err = ", err)
#     }
#   }
#   if (ii == max_iter) return(-Inf)
# }

