# ---- common functions ---- 

MaxRelErr <- function(lambdaNew, lambdaOld) {
  # lambdaSum <- abs(lambdaNew) + abs(lambdaOld)
  # inds <- (lambdaSum < 1e-10) # indices of entires that are very close to 0
  # relErr <- rep(NA,length(lambdaNew))
  # # those entires that are very close to 0 use abs err
  # relErr[inds] <- abs(lambdaNew - lambdaOld)[inds]
  # # other entires use relative err
  # relErr[!inds] <- abs((lambdaNew - lambdaOld) / (lambdaNew + lambdaOld))[!inds]
  
  # TODO: added for numerical stability, what is a good tolerance ?
  if (max(abs(lambdaNew - lambdaOld)) < 1e-10) return(0)
  relErr <- abs((lambdaNew - lambdaOld) / (lambdaNew + lambdaOld))
  return(max(relErr))
}

# adjusted EM by chen-et-al2008
adjG_R <- function(G, an) {
  n <- nrow(G)
  gbar <- 1/n*colSums(G)
  if (missing(an)) an <- max(1,0.5*log(n))
  gadd <- -an*gbar
  return(rbind(G,gadd))
}

# ---- non-censoring EL ----

log.star <- function(x, n) {
  cond <- x >= 1/n
  ans <- rep(NaN, length(x))
  ans[cond] <- log(x[cond])
  ans[!cond] <- -1/2 * n^2 * x[!cond]^2 + 2*n*x[!cond] - 3/2 - log(n)
  ans
}

log.star1 <- function(x, n) {
  cond <- x >= 1/n
  ans <- rep(NaN,length(x))
  ans[cond] <- 1/(x[cond])
  ans[!cond] <- -n^2*x[!cond] + 2*n
  return(ans)
}

log.star2 <- function(x, n) {
  cond <- x >= 1/n
  ans <- rep(NaN,length(x))
  ans[cond] <- -1/(x[cond]^2)
  ans[!cond] <- -n^2
  return(ans)
}

# Note: G is nObs x nEqs
Qfun <- function(lambda, G) {
  G <- t(G)
  N <- ncol(G) # nObs
  sum(apply(G, 2, function(gg) log.star(x = 1 - sum(lambda*gg), n = N)))
}

# R implementation of lambdaNR
# Input: G is nObs x nEqs matrix
lambdaNR_R <- function(G, max_iter=100, rel_tol=1e-7, verbose = FALSE, 
                       lambdaOld = NULL) {
  G <- t(G)
  nObs <- ncol(G)
  nEqs <- nrow(G)
  if (is.null(lambdaOld)) {
    # lambdaOld <- rnorm(nEqs) # TODO: initialize to random value?
    lambdaOld <- rep(0,nEqs)
  }
  lambdaNew <- lambdaOld
  nIter <- 0
  for (ii in 1:max_iter) {
    nIter <- ii
    Glambda <- t(lambdaOld) %*% G
    Glambda <- 1 - Glambda
    Q2 <- matrix(0, nrow=nEqs, ncol=nEqs)
    rho <- rep(NaN, nObs)
    for(jj in 1:nObs) {
      rho[jj] <- log.star1(Glambda[jj],nObs)
      Q2 <- Q2 + log.star2(Glambda[jj],nObs) * (G[,jj] %*% t(G[,jj]))
    }
    Q1 <- -G %*% rho
    lambdaNew <- lambdaOld - solve(Q2,Q1) # TODO: step size?
    maxErr <- MaxRelErr(lambdaNew, lambdaOld)
    # if (verbose && (nIter %% 5 == 0)){
    if (verbose){
      message("nIter = ", nIter)
      message("err = ", maxErr)
    }
    if (maxErr < rel_tol) break;
    lambdaOld <- lambdaNew
  }
  notconv <- (ii == max_iter && maxErr > rel_tol)
  if(notconv) lambdaNew <- rep(NaN, nEqs)
  # c(lambdaNew) to make sure lambda is a vector
  output <- list(lambda=c(lambdaNew), convergence=!notconv)
  return(output)
}

# G is nObs x nEqs matrix
omega.hat.NC_R <- function(G, adjust = FALSE, 
                           max_iter = 100, rel_tol = 1e-07, verbose = FALSE) {
  lambdaOut <- lambdaNR_R(G = G, max_iter, rel_tol, verbose)
  # nIter <- lambdaout$nIter
  # maxErr <- lambdaout$maxErr
  # notconv <- nIter == max_iter && maxErr > rel_tol # 1 if not converged
  conv <- lambdaOut$convergence # 1 if converged
  if (!conv) {
    nObs <- nrow(G)
    omegahat <- rep(NaN,nObs)
  }
  else {
    lambdahat <- lambdaOut$lambda
    omegahat <- c(1/(1-t(lambdahat) %*% t(G)) / sum(1/(1-t(lambdahat) %*% t(G))))
  }
  # returns a vector of omegahat
  return(omegahat)
  # return(list(omegas=omegahat, convergence=conv))
}

#---- censoring EL ----

# log.sharp and related are unique to censoring
# Note: x and q must be of the same length
log.sharp <- function(x, q) {
  cond <- x >= q
  ans <- rep(NaN,length(x))
  ans[cond] <- log(x[cond])
  ans[!cond] <- -1/(2*q[!cond]^2)*x[!cond]^2 + 2/q[!cond]*x[!cond] - 3/2 + log(q[!cond])
  return(ans)
}

log.sharp1 <- function(x, q) {
  cond <- x >= q
  ans <- rep(NaN,length(x))
  ans[cond] <- 1/(x[cond])
  ans[!cond] <- -1/(q[!cond]^2)*x[!cond] + 2/q[!cond]
  return(ans)
}

log.sharp2 <- function(x, q) {
  cond <- x >= q
  ans <- rep(NaN,length(x))
  ans[cond] <- -1/(x[cond]^2)
  ans[!cond] <- -1/q[!cond]^2
  return(ans)
}

# for mle.check
# Input: G is nObs x nEqs
QfunCens <- function(lambda, G, weights) {
  G <- t(G)
  weights_sum <- sum(weights) # this should be nObs tho..
  G_list <- split(G, rep(1:ncol(G), each = nrow(G)))
  ls <- mapply(function(gg,qq) log.sharp(x = weights_sum + sum(lambda*gg), qq),
               G_list, weights)
  # as.numeric otherwise it's a 1-1 matrix
  return(as.numeric(weights %*% ls))
}

# R implementation of lambdaNRC
# G is nObs x nEqs matrix 
# lambdaOld passed in is for the EM algorithm `omega.hat.EM_R` to work properly
lambdaNRC_R <- function(G, weights, max_iter = 100, rel_tol = 1e-7, verbose = FALSE,
                        lambdaOld = NULL) {
  G <- t(G)
  nObs <- ncol(G)
  nEqs <- nrow(G)
  if (is.null(lambdaOld)) lambdaOld <- rep(0,nEqs)
  lambdaNew <- lambdaOld
  nIter <- 0
  # newton-raphson loop
  for (ii in 1:max_iter) {
    nIter <- ii
    # Q1 and Q2
    Glambda <- t(lambdaOld) %*% G
    Glambda <- sum(weights) + Glambda
    # Q1 <- rep(0,nEqs)
    rho <- rep(NaN,nObs)
    Q2 <- matrix(rep(0,nEqs*nEqs), nEqs, nEqs)
    for (jj in 1:nObs) {
      rho[jj] <- log.sharp1(Glambda[jj], weights[jj])
      # to avoid numerical problem?
      if (weights[jj] > rel_tol*nObs) {
        Q2 <- Q2 + weights[jj]*log.sharp2(Glambda[jj], weights[jj])*(G[,jj] %*% t(G[,jj]))
      }
    }
    weights.sav <- weights
    weights.sav[weights.sav < rel_tol*nObs] <- 0
    rho[is.infinite(rho) & !weights.sav] <- 0 
    # Q1 <- G %*% (rho * weights)
    Q1 <- G %*% (rho * weights.sav)
    # print("Q2invQ1 = ")
    # print(t(solve(Q2,Q1)))
    lambdaNew <- lambdaOld - solve(Q2,Q1)
    maxErr <- MaxRelErr(lambdaNew, lambdaOld) # maximum relative error
    # message("maxErr = ", maxErr)
    if (maxErr < rel_tol) {
      break;
    }
    lambdaOld <- lambdaNew # complete cycle
    if (verbose && (nIter %% 10 == 0)){
      message("nIter = ", nIter)
      message("err = ", maxErr)
    }
  }
  notconv <- (ii == max_iter && maxErr > rel_tol)
  if (notconv) lambdaNew <- rep(NaN, nEqs)
  output<- list(lambda = c(lambdaNew), convergence=!notconv)
  return(output)
}

evalPsos_R <- function(ii, epsOrd, omegas) {
  nObs <- length(omegas)
  psos <- 0
  for (jj in nObs:1) {
    kk <- epsOrd[jj]
    psos <- psos + omegas[kk] # sum over the omegas of eps <= than ii-th eps
    if (kk == ii) break
  }
  return(psos)
}

# calculate the weights for the weighted (censored) EL problem 
evalWeights_R <- function(deltas, omegas, epsilons) {
  nObs <- length(omegas)
  epsOrd <- order(epsilons)
  psots <- rep(0,nObs)
  for (ii in 1:nObs) {
    for (jj in 1:nObs) {
      kk <- epsOrd[jj]
      if (deltas[kk] == 0) {
        psots[ii] <- psots[ii] + omegas[ii]/evalPsos_R(kk, epsOrd, omegas)
      }
      if (kk == ii) break
    }
  }
  # the weights have the order as the original sample 
  weights <- deltas + psots
  return(weights)
}

evalEpsilons_R <- function(y,X,beta) {
  eps <- y - X %*% beta
  return(eps)
}

evalEpsilonsLS_R <- function(y,X,Z,beta,gamma,sig2) {
  eps <- (y - X %*% beta)*exp(-Z %*% gamma)/sqrt(sig2)
  return(eps)
}

# G is nObs x nEqs matrix 
omega.hat.EM_R <- function(G, deltas, epsilons, adjust = FALSE, 
                           max_iter = 100, rel_tol = 1e-7, verbose=FALSE) {
  n <- nrow(G)
  m <- ncol(G)
  err <- Inf
  # lambdaOld <- rep(0,m)
  nIter <- 0
  # initialize omegas with uncensored solution 
  omegas <- omega.hat.NC_R(G, adjust, max_iter, rel_tol, verbose)
  if (any(is.nan(omegas))) {
    message("Initial omegas are nans.")
    return(rep(NaN,length(deltas)))
  }
  # if (adjust) omegas <- omegas[1:(n-1)]
  if (adjust) {
    epsilons <- c(epsilons,0)
    deltas <- c(deltas,0)
  }
  # omegasOld <- omegas
  logelOld <- logEL_R(omegas,epsilons,deltas)
  for (ii in 1:max_iter) {
    nIter <- ii
    # E step: calculating weights
    weights <- evalWeights_R(deltas, omegas, epsilons)
    # if (adjust) weights <- c(weights,0)
    # M step:
    # lambdaOut <- lambdaNRC_R(G, weights, max_iter, rel_tol, verbose, lambdaOld)
    lambdaOut <- lambdaNRC_R(G, weights, max_iter, rel_tol, verbose)
    # TODO: what if not converged ?? use a random weights and continue ?
    if (!lambdaOut$convergence) {
      # message("lambdaNRC did not converge in EM")
      return(rep(NaN,n))
    }
    lambdaNew <- lambdaOut$lambda
    qlg <- c(sum(weights) + lambdaNew %*% t(G))
    omegas <- weights/qlg
    omegas <- omegas/sum(omegas)
    if (any(omegas < -rel_tol)) message("omega.hat.EM_R: negative omegas.")
    omegas <- abs(omegas)
    # err <- MaxRelErr(lambdaNew,lambdaOld)
    # err <- MaxRelErr(omegas,omegasOld)
    # if (adjust) omegas <- omegas[1:(n-1)]
    logel <- logEL_R(omegas,epsilons,deltas)
    err <- MaxRelErr(logel,logelOld)
    if (verbose && nIter %% 20 == 0) {
      message("nIter = ", nIter)
      message("err = ", err)
    }
    if (err < rel_tol) break
    # lambdaOld <- lambdaNew
    # omegasOld <- omegas
    logelOld <- logel
  }
  notconv <- (nIter == max_iter && err > rel_tol) # TRUE if not converged 
  if (notconv) {
    omegas = rep(NaN,n)
  }
  return(omegas)
}

# accelerate EM wang-et-al08
vecinv_R <- function(x) {
  return(x/sum(x*x))
}


# omega.hat.EM_R <- function(G, deltas,
#                            max_iter = 100, rel_tol = 1e-7, verbose=FALSE) {
#   n <- nrow(G)
#   m <- ncol(G)
#   err <- Inf
#   # lambdaOld <- rep(0,m)
#   nIter <- 0
#   # initialize omegas with uncensored solution 
#   omegas <- omega.hat.NC_R(G, max_iter, rel_tol, verbose)
#   if (any(is.nan(omegas))) {
#     message("Initial omegas are nans.")
#     return(rep(NaN,length(deltas)))
#   }
#   # omegasOld <- omegas
#   logelOld <- logEL_R(omegas,epsilons,deltas)
#   if (any(is.nan(omegas))) return(rep(NaN,n))
#   for (ii in 1:max_iter) {
#     nIter <- ii
#     # E step: calculating weights
#     weights <- evalWeights_R(deltas, omegas, epsilons)
#     # M step:
#     # lambdaOut <- lambdaNRC_R(G, weights, max_iter, rel_tol, verbose, lambdaOld)
#     lambdaOut <- lambdaNRC_R(G, weights, max_iter, rel_tol, verbose)
#     # TODO: what if not converged ?? use a random weights and continue ?
#     if (!lambdaOut$convergence) {
#       # message("lambdaNRC did not converge in EM")
#       return(rep(NaN,n))
#     }
#     lambdaNew <- lambdaOut$lambda
#     qlg <- c(sum(weights) + lambdaNew %*% t(G))
#     omegas <- weights/qlg
#     omegas <- omegas/sum(omegas)
#     # err <- MaxRelErr(lambdaNew,lambdaOld)
#     # err <- MaxRelErr(omegas,omegasOld)
#     logel <- logEL_R(omegas,epsilons,deltas)
#     err <- MaxRelErr(logel,logelOld)
#     if (verbose && nIter %% 20 == 0) {
#       message("nIter = ", nIter)
#       message("err = ", err)
#     }
#     if (err < rel_tol) break
#     # lambdaOld <- lambdaNew
#     # omegasOld <- omegas
#     logelOld <- logel
#   }
#   notconv <- (nIter == max_iter && err > rel_tol) # TRUE if not converged 
#   if (notconv) {
#     omegas = rep(NaN,n)
#   }
#   return(omegas)
# }

# ---- wrapper functions ----

# wrapper function for censor and non-censor omega.hat
omega.hat_R <- function(G, deltas, epsilons, adjust=FALSE,
                        max_iter = 100, rel_tol = 1e-7, verbose = FALSE) {
  if (missing(deltas) && missing(epsilons)) {
    omegas <- omega.hat.NC_R(G, adjust, max_iter=max_iter, rel_tol=rel_tol, verbose=verbose)
  }
  else {
    omegas <- omega.hat.EM_R(G, deltas, epsilons, adjust, 
                             max_iter=max_iter, rel_tol=rel_tol, verbose=verbose)
  }
  return(c(omegas))
}

# G is nObs x nEqs matrix
logEL_R <- function(omegas, epsilons, deltas, adjust=FALSE) {
  if (any(is.nan(omegas))) return(-Inf)
  if (missing(deltas)) {
    sum(log(omegas))
  }
  else {
    if (adjust) {
      epsilons <- c(epsilons,0)
      deltas <- c(deltas,0)
    }
    epsOrd <- order(epsilons) # ascending order of epsilons
    n <- length(omegas)
    psos <- rep(0,n)
    for (ii in 1:n) {
      psos[ii] <- evalPsos_R(ii, epsOrd, omegas) 
    }
    # numerical stability: watch out for extremely small negative values
    omegas[abs(omegas) < 1e-10/length(omegas)] <- 1e-10
    return(sum(deltas*log(omegas)+(1-deltas)*log(psos)))
  }
}

# # G is nObs x nEqs matrix
# logEL_R <- function(G, deltas, epsilons, max_iter = 100, rel_tol = 1e-7) {
#   # non-censored case:
#   if (missing(deltas) && missing(epsilons)) {
#     omegas <- omega.hat.NC_R(G, max_iter, rel_tol)
#     if (any(is.nan(omegas))) return(-Inf)
#     else return(sum(log(omegas)))
#   }
#   # censored case:
#   else {
#     if (length(deltas) != length(epsilons)) {
#       stop("deltas and epsilons have inconsistent dimensions.")
#     }
#     if (nrow(G) != length(deltas)) {
#       stop("deltas and G have inconsistent dimensions.")
#     }
#     omegas <- omega.hat.EM_R(G, deltas, epsilons, max_iter, rel_tol)
#     if (any(is.nan(omegas))) return(-Inf)
#     else {
#       epsOrd <- order(epsilons) # ascending order of epsilons
#       n <- length(omegas)
#       psos <- rep(0,n)
#       for (ii in 1:n) {
#         psos[ii] <- evalPsos_R(ii, epsOrd, omegas) 
#       }
#       return(sum(deltas*log(omegas)+(1-deltas)*log(psos)))
#     }
#   }
# }
