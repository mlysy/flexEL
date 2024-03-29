# ---- helper functions ----

# wrappers to R versions having the same interface as C++

lambda_nr_R <- function(G, max_iter = 100, rel_tol = 1e-7,
                        supp_adj = FALSE, supp_adj_a, weights = NULL) {
  if(supp_adj) {
    G <- adjG_R(G, an = supp_adj_a)
  }
  if(is.null(weights)) {
    ans <- lambdaNR_R(G, max_iter = max_iter, rel_tol = rel_tol)
  } else {
    if(supp_adj) weights <- c(weights, 1)
    weights <- weights/sum(weights)
    ans <- lambdaNRC_R(G, weights = weights,
                       max_iter = max_iter, rel_tol = rel_tol)
  }
  ans
}

omega_hat_R <- function(G, max_iter = 100, rel_tol = 1e-7,
                        supp_adj = FALSE, supp_adj_a, weights = NULL) {
  if(supp_adj) {
    G <- adjG_R(G, an = supp_adj_a)
  }
  if(is.null(weights)) {
    ans <- omega_hat_NC_R(G, max_iter = max_iter, rel_tol = rel_tol)
  } else {
    if(supp_adj) weights <- c(weights, 1)
    weights <- weights/sum(weights)
    ans <- weighted_omega_hat_R(G, weights = weights,
                                max_iter = max_iter, rel_tol = rel_tol)
  }
  ans
}

logel_R <- function(G, max_iter = 100, rel_tol = 1e-7,
                    supp_adj = FALSE, supp_adj_a, weights = NULL) {
  if(supp_adj) {
    G <- adjG_R(G, an = supp_adj_a)
  }
  if(is.null(weights)) {
    omega_hat <- omega_hat_NC_R(G, max_iter = max_iter, rel_tol = rel_tol)
    lel <- if(anyNA(omega_hat)) -Inf else sum(log(omega_hat))
  } else {
    if(supp_adj) weights <- c(weights, 1)
    lel <- weighted_logEL_G_R(G, weights = weights,
                              max_iter = max_iter, rel_tol = rel_tol)
    lel <- as.numeric(lel)
  }
  lel
}

logel_grad_R <- function(G, max_iter = 100, rel_tol = 1e-7,
                         supp_adj = FALSE, supp_adj_a, weights = NULL) {
  n <- nrow(G)
  p <- ncol(G)
  fun <- function(G) {
    logel_R(G, weights = weights,
            max_iter = max_iter, rel_tol = rel_tol,
            supp_adj = supp_adj, supp_adj_a = supp_adj_a)
  }
  dldG <- tryCatch(numDeriv::grad(fun, G),
                   error = function(e) rep(NaN, n*p))
  matrix(dldG, nrow = n, ncol = p)
}


# maximum relative error
max_rel_err <- function(lambdaNew, lambdaOld) {
  # Note: for stability (using the same tolerance as in C++ implementation)
  ## if (max(abs(lambdaNew - lambdaOld)) < 1e-10) return(0)
  ## relErr <- abs((lambdaNew - lambdaOld) / (lambdaNew + lambdaOld))
  ## return(max(relErr))
  max(abs(lambdaNew - lambdaOld)/(.1 + abs(lambdaNew + lambdaOld)))
}

# adjusted EM by chen-et-al2008
adjG_R <- function(G, an) {
  n <- nrow(G)
  gbar <- 1/n*colSums(G)
  if (missing(an)) an <- max(1,0.5*log(n))
  gadd <- -an*gbar
  return(unname(rbind(G,gadd)))
}

# ---- EL (no censor) ----

log_star <- function(x, n) {
  cond <- x >= 1/n
  ans <- rep(NaN, length(x))
  ans[cond] <- log(x[cond])
  ans[!cond] <- -1/2 * n^2 * x[!cond]^2 + 2*n*x[!cond] - 3/2 - log(n)
  ans
}

# 1st derivative of log_star
log_star1 <- function(x, n) {
  cond <- x >= 1/n
  ans <- rep(NaN,length(x))
  ans[cond] <- 1/(x[cond])
  ans[!cond] <- -n^2*x[!cond] + 2*n
  return(ans)
}

# 2nd derivative of log_star
log_star2 <- function(x, n) {
  cond <- x >= 1/n
  ans <- rep(NaN,length(x))
  ans[cond] <- -1/(x[cond]^2)
  ans[!cond] <- -n^2
  return(ans)
}

# G is nObs x nEqs
Qfun <- function(lambda, G) {
  G <- t(G)
  N <- ncol(G) # nObs
  sum(apply(G, 2, function(gg) log_star(x = 1 - sum(lambda*gg), n = N)))
}

# R implementation of lambdaNR
# G is nObs x nEqs
lambdaNR_R <- function(G, max_iter=100, rel_tol=1e-7, verbose = FALSE,
                       lambdaOld = NULL) {
  G <- t(G)
  nObs <- ncol(G)
  nEqs <- nrow(G)
  if (is.null(lambdaOld)) {
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
      rho[jj] <- log_star1(Glambda[jj],nObs)
      Q2 <- Q2 + log_star2(Glambda[jj],nObs) * (G[,jj] %*% t(G[,jj]))
    }
    Q1 <- -G %*% rho
    lambdaNew <- lambdaOld - solve(Q2,Q1) # TODO: step size?
    maxErr <- max_rel_err(lambdaNew, lambdaOld)
    # if (verbose && (nIter %% 5 == 0)){
    if (verbose) {
      message("nIter = ", nIter)
      message("err = ", maxErr)
    }
    if (maxErr < rel_tol) break;
    lambdaOld <- lambdaNew
  }
  notconv <- (ii == max_iter && maxErr > rel_tol)
  if(notconv) lambdaNew <- rep(NaN, nEqs)
  # Note: c(lambdaNew) to make sure lambda is a vector
  output <- list(lambda=c(lambdaNew), convergence=!notconv)
  return(output)
}

# G is nObs x nEqs
omega_hat_NC_R <- function(G, max_iter = 100, rel_tol = 1e-7, verbose = FALSE) {
  lambdaOut <- lambdaNR_R(G = G, max_iter, rel_tol, verbose)
  conv <- lambdaOut$convergence # 1 if converged
  # message("lambda = ")
  # print(lambdaOut$lambda)
  if (!conv) {
    nObs <- nrow(G)
    omegahat <- rep(NaN, nObs)
  }
  else {
    lambdahat <- lambdaOut$lambda
    omegahat <- c(1/(1-t(lambdahat) %*% t(G)) / sum(1/(1-t(lambdahat) %*% t(G))))
  }
  # returns a vector of omegahat
  return(omegahat)
}

weighted_omega_hat_R <- function(G, weights, adjust = FALSE,
                                 max_iter = 100, rel_tol = 1e-7, abs_tol = 1e-3,
                                 verbose = FALSE) {
  lambdaOut <- lambdaNRC_R(G = G, weights = weights, max_iter, rel_tol, verbose)
  conv <- lambdaOut$convergence # 1 if converged
  # message("lambda = ")
  # print(lambdaOut$lambda)
  if (!conv) {
    nObs <- nrow(G)
    omegahat <- rep(NaN, nObs)
  }
  else {
    lambdahat <- lambdaOut$lambda
    omegahat <- c(weights/(1-lambdahat %*% t(G)))
    # omegahat <- c(1/(1-t(lambdahat) %*% t(G)) / sum(1/(1-t(lambdahat) %*% t(G))))
  }
  # returns a vector of omegahat
  return(omegahat)
}

# ---- EL (with censor) ----

# Note: log_sharp and related are unique to censoring
log_sharp <- function(x, q) {
  if (length(x) != length(q)) stop("x and q must be of the same length.")
  cond <- x >= q
  ans <- rep(NaN,length(x))
  ans[cond] <- log(x[cond])
  ans[!cond] <- -1/(2*q[!cond]^2)*x[!cond]^2 + 2/q[!cond]*x[!cond] - 3/2 + log(q[!cond])
  return(ans)
}

# 1st derivative of log_sharp
log_sharp1 <- function(x, q) {
  cond <- x >= q
  ans <- rep(NaN,length(x))
  ans[cond] <- 1/(x[cond])
  ans[!cond] <- -1/(q[!cond]^2)*x[!cond] + 2/q[!cond]
  return(ans)
}

# 2nd derivative of log_sharp
log_sharp2 <- function(x, q) {
  cond <- x >= q
  ans <- rep(NaN,length(x))
  ans[cond] <- -1/(x[cond]^2)
  ans[!cond] <- -1/q[!cond]^2
  return(ans)
}

# G is nObs x nEqs
QfunCens <- function(lambda, G, weights) {
  G <- t(G)
  weights_sum <- sum(weights) # this should be nObs tho..
  G_list <- split(G, rep(1:ncol(G), each = nrow(G)))
  ls <- mapply(function(gg,qq) log_sharp(x = weights_sum + sum(lambda*gg), qq),
               G_list, weights)
  # as.numeric otherwise it's a 1-1 matrix
  return(as.numeric(weights %*% ls))
}

# R implementation of lambdaNRC (weighted lambdaNR)
# G is nObs x nEqs
# weights should be normalized (sum to 1)
# lambdaOld passed in is for the EM algorithm `_R` to work properly
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
    Glambda <- sum(weights) - Glambda
    # Q1 <- rep(0,nEqs)
    rho <- rep(NaN,nObs)
    Q2 <- matrix(rep(0,nEqs*nEqs), nEqs, nEqs)
    for (jj in 1:nObs) {
      rho[jj] <- log_sharp1(Glambda[jj], weights[jj])
      # to avoid numerical problem? (begin)
      # if (weights[jj] > rel_tol*nObs) {
      #   Q2 <- Q2 - weights[jj]*log_sharp2(Glambda[jj], weights[jj])*(G[,jj] %*% t(G[,jj]))
      # }
      # (end)
      Q2 <- Q2 - weights[jj]*log_sharp2(Glambda[jj], weights[jj])*(G[,jj] %*% t(G[,jj]))
    }
    # to avoid numerical problem? (begin)
    # weights_sav <- weights
    # weights_sav[weights_sav < rel_tol*nObs] <- 0
    # rho[is.infinite(rho) & !weights_sav] <- 0
    # Q1 <- G %*% (rho * weights_sav)
    # (end)
    Q1 <- G %*% (rho * weights)
    lambdaNew <- lambdaOld - solve(Q2,Q1)
    maxErr <- max_rel_err(lambdaNew, lambdaOld) # maximum relative error
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

# evaluate Partial Sum of OmegaS
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
  nObs <- length(deltas)
  epsOrd <- order(epsilons)
  # message("epsOrd = ")
  # print(epsOrd-1)
  psots <- rep(0,nObs)
  for (ii in 1:nObs) {
    for (jj in 1:nObs) {
      kk <- epsOrd[jj]
      if (deltas[kk] == 0) {
        # print(evalPsos_R(kk, epsOrd, omegas))
        psots[ii] <- psots[ii] + omegas[ii]/evalPsos_R(kk, epsOrd, omegas)
      }
      if (kk == ii) break
    }
  }
  # the weights have the order as the original sample
  weights <- deltas + psots
  return(weights)
}

# G is nObs x nEqs
# Note: if adjust G (using support adjustment), pass the adjusted G to this function and set adjust to TRUE
omega_hat_EM_R <- function(G, deltas, epsilons, adjust = FALSE,
                           max_iter = 100, rel_tol = 1e-7, abs_tol = 1e-3, verbose=FALSE,
                           dbg = FALSE) {
  n <- nrow(G)
  m <- ncol(G)
  err <- Inf
  # lambdaNew <- rep(0,m)
  nIter <- 0
  # # initialize omegas with uncensored solution
  # omegas <- omega_hat_NC_R(G, max_iter, rel_tol, verbose=FALSE)
  omegas <- rep(1/n, n)
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
  logelOld <- logEL_R(omegas, epsilons, deltas)
  # message("logelOld = ", logelOld)

  # weights <- evalWeights_R(deltas, omegas, epsilons)
  # weights <- weights/sum(weights)
  # print(weights/sum(weights))
  # message("logelOld2 = ", sum(weights*log(omegas)))

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
    weights <- evalWeights_R(deltas, omegas, epsilons)
    weights <- weights/sum(weights)
    # print(weights)
    # if (ii == 1) {
    #   print(weights)
    # }
    # M step:
    lambdaOut <- lambdaNRC_R(G, weights, max_iter, rel_tol, verbose=FALSE)
    # TODO: what if not converged ?? use a random weights and continue ?
    if (!lambdaOut$convergence) {
      message("lambdaNRC did not converge in EM")
      if (dbg) {
        return(list(conv=FALSE, omegas=rep(NaN,n), logels=logels, artome=artome))
      }
      else {
        return(list(conv=FALSE, omegas=rep(NaN,n)))
      }
    }
    lambdaNew <- lambdaOut$lambda
    qlg <- c(sum(weights) - lambdaNew %*% t(G))
    omegas <- weights/qlg
    # omegas <- omegas/sum(omegas)
    if (any(omegas < -rel_tol)) {
      message("omega_hat_EM_R: negative omegas.")
      message(print(omegas))
    }
    omegas <- abs(omegas)
    omegas <- omegas/sum(omegas)
    # err <- max_rel_err(lambdaNew,lambdaOld)
    # err <- max_rel_err(omegas,omegasOld)
    logel <- logEL_R(omegas,epsilons,deltas)
    # message("logel = ", logel)

    # for debug
    if (dbg) {
      omegamMat[,ii+1] <- omegas
      logels <- c(logels,logel)
      artome <- c(artome,weights[n])
    }

    # err <- max_rel_err(logel,logelOld)
    err <- abs(logel-logelOld)
    # if (verbose && nIter %% 20 == 0) {
    if (verbose) {
      message("nIter = ", nIter)
      # message("err = ", err)
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

# ---- EL (with continuity correction) ----

# smoothed indicator function 1(x <= 0)
ind_smooth_R <- function(x, s=10) {
  return(1/(1+exp(s*x)))
}

# 1st derivative of smoothed indicator function 1(x <= 0)
ind1_smooth_R <- function(x, s=10) {
  return(-s*exp(s*x)/(1+exp(s*x))^2)
}

# smoothed partial sum of omegas (vectorized version)
evalPsos_smooth_R <- function(ii, omegas, epsilons, s = 10, support = FALSE) {
  n <- length(epsilons)
  if (support == TRUE  && ii == n) {
    psos <- sum(omegas[1:(n-1)])+0.5*omegas[n]
  }
  else psos <- sum(ind_smooth_R(epsilons[ii]-epsilons,s)*omegas)
  return(psos)
}

# smoothed censored logEL
logEL_smooth_R <- function(omegas, epsilons, deltas, s = 10, adjust=FALSE) {
  if (adjust) {
    epsilons <- c(epsilons,-Inf)
    deltas <- c(deltas,0)
  }
  epsOrd <- order(epsilons) # ascending order of epsilons
  n <- length(omegas)
  psos <- rep(0,n)
  for (ii in 1:n) {
    psos[ii] <- evalPsos_smooth_R(ii,omegas,epsilons,s,support = adjust)
  }
  # numerical stability: watch out for extremely small negative values
  omegas[abs(omegas) < 1e-10/length(omegas)] <- 1e-10

  weights <- evalWeights_smooth_R(deltas, omegas, epsilons, s, adjust)
  return(sum(weights*log(omegas)))
  # return(sum(deltas*log(omegas)+(1-deltas)*log(psos)))
}

# smoothed evalWeight function with sigmoid function (more vectorized version)
evalWeights_smooth_R <- function(deltas, omegas, epsilons, s=10, support = FALSE) {
  n <- length(epsilons)
  epsOrd <- order(epsilons)
  psots <- rep(0,n)
  psoss <- rep(NA,n)
  for (jj in 1:n) {
    psoss[jj] <- evalPsos_smooth_R(jj,omegas,epsilons,s,support)
  }
  # print(psoss)
  if (support) {
    for (ii in 1:n) {
      temp <- ind_smooth_R(epsilons-epsilons[ii],s)
      if (ii==n) temp[n] <- ind_smooth_R(0,s)
      psots[ii] <- sum((1-deltas)*temp*omegas[ii]/psoss)
      # message(psots[ii])
    }
  }
  else {
    for (ii in 1:n) {
      psots[ii] <- sum((1-deltas)*ind_smooth_R(epsilons-epsilons[ii],s)*omegas[ii]/psoss)
    }
  }
  weights <- deltas + psots
  return(weights)
}

# using smoothed objective function in EM
omega_hat_EM_smooth_R <- function(G, deltas, epsilons, s=10, adjust = FALSE,
                                  max_iter = 100, rel_tol = 1e-7, abs_tol = 1e-3,
                                  verbose=FALSE, dbg = FALSE) {
  n <- nrow(G)
  m <- ncol(G)
  err <- Inf
  nIter <- 0
  # initialize omegas with uncensored solution
  # omegas <- omega_hat_NC_R(G, max_iter, rel_tol, verbose=FALSE)
  omegas <- rep(1/n, n)
  if (any(is.nan(omegas))) {
    message("Initial omegas are nans.")
    # return(rep(NaN,length(deltas)))
    return(list(conv=FALSE, omegas=rep(NaN,length(deltas))))
  }
  # if (adjust) {
  #   epsilons <- c(epsilons,-Inf)
  #   deltas <- c(deltas,0)
  # }
  logelOld <- logEL_smooth_R(omegas,epsilons,deltas,s,adjust = adjust)

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
    if (adjust) {
      weights <- evalWeights_smooth_R(c(deltas,0), omegas, c(epsilons, -Inf), s, support = adjust)
    } else {
      weights <- evalWeights_smooth_R(deltas, omegas, epsilons, s, support = adjust)
    }
    # message(paste0(weights, collapse = " "))
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
    qlg <- c(sum(weights) - lambdaNew %*% t(G))
    omegas <- weights/qlg
    # omegas <- omegas/sum(omegas)
    if (any(omegas < -rel_tol)) message("omega_hat_EM_R: negative omegas.")
    omegas <- abs(omegas)
    omegas <- omegas/sum(omegas)
    logel <- logEL_smooth_R(omegas,epsilons,deltas,s,adjust = adjust)

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

# ---- wrapper functions ----

# wrapper function for censor and non-censor omega_hat
omega_hat_cens_R <- function(G, deltas, epsilons, adjust = FALSE,
                        max_iter = 100, rel_tol = 1e-7, abs_tol = 1e-3, verbose = FALSE) {
  if (missing(deltas) && missing(epsilons)) {
    omegas <- omega_hat_NC_R(G, max_iter=max_iter, rel_tol=rel_tol, verbose=verbose)
  }
  else {
    omegas <- omega_hat_EM_R(G, deltas, epsilons, adjust,
                             max_iter=max_iter, rel_tol=rel_tol, abs_tol = abs_tol, verbose=verbose)
  }
  return(c(omegas))
}

# G is nObs x nEqs
logEL_R <- function(omegas, epsilons, deltas, adjust=FALSE) {
  if (any(is.nan(omegas))) return(-Inf)
  if (missing(deltas)) {
    return(sum(log(omegas)))
  }
  else {
    if (adjust) {
      epsilons <- c(epsilons,-Inf)
      deltas <- c(deltas,0)
    }
    # if (adjust) omegas <- omegas[1:length(epsilons)]
    epsOrd <- order(epsilons) # ascending order of epsilons
    n <- length(omegas)
    psos <- rep(0,n)
    for (ii in 1:n) {
      psos[ii] <- evalPsos_R(ii, epsOrd, omegas)
    }
    # numerical stability: watch out for extremely small negative values
    omegas[abs(omegas) < 1e-10/length(omegas)] <- 1e-10

    weights <- evalWeights_R(deltas, omegas, epsilons)
    # print(weights)
    # message("new weights logel: ", sum(weights*log(omegas)))
    return(sum(weights*log(omegas)))
    # return(sum(deltas*log(omegas)+(1-deltas)*log(psos)))
  }
}

logEL_G_R <- function(G, max_iter = 100, rel_tol = 1e-7) {
  # lambda <- flexEL::lambdaNR(G = G, max_iter = max_iter, rel_tol = rel_tol, support = FALSE)
  lambda_lst <- lambdaNR_R(G, max_iter, rel_tol)
  if (!lambda_lst$convergence) return(NA)
  else {
    lambda <- lambda_lst$lambda
    omega <- 1/nrow(G) * 1/(1 - c(G %*% lambda))
    neglogel <- sum(log(omega))
    dldG <- length(omega) * t(lambda %*% t(omega))
    attr(neglogel, "gradient") <- dldG
    return(neglogel)
  }
}

weighted_logEL_G_R <- function(G, weights, max_iter = 100, rel_tol = 1e-7) {
  sum_weights <- sum(weights)
  norm_weights <- weights/sum_weights
  lambda_lst <- lambdaNRC_R(G, norm_weights, max_iter, rel_tol)
  if (!lambda_lst$convergence) return(-Inf)
  else {
    lambda <- lambda_lst$lambda
    omega <- c(norm_weights/(1 - c(G %*% lambda)))
    logel <- sum(norm_weights * log(omega)) * sum_weights
    dldG <- t(lambda %*% t(omega)) * sum_weights
    attr(logel, "gradient") <- dldG
    return(logel)
  }
}

logEL_adjG_R <- function(G, max_iter = 100, rel_tol = 1e-7,
                         an = max(1, 0.5*log(nrow(G)))) {
  n <- nrow(G)
  # lambda <- flexEL::lambdaNR(G = G, max_iter = max_iter, rel_tol = rel_tol, support = TRUE)
  # omega <- flexEL:::omega_hat(G = G, max_iter = max_iter, rel_tol = rel_tol, support = TRUE)
  lambda <- lambdaNR_R(adjG_R(G, an), max_iter = max_iter, rel_tol = rel_tol)$lambda
  G_adj <- adjG_R(G, an)
  omega <- c(1/(1-t(lambda) %*% t(G_adj)) / sum(1/(1-t(lambda) %*% t(G_adj))))
  neglogel <- sum(log(omega))
  # an <- max(1, 0.5*log(n))
  dldGadj <- (n+1) * t(lambda %*% t(omega[1:n])) - (n+1) * omega[n+1]*an/n * rep(1, n) %*% t(lambda)
  attr(neglogel, "gradient") <- dldGadj
  return(neglogel)
}

