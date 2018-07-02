# ---- mean regression ----

# location model
# mr.evalG is the same for non-censoring and censoring for mean regression
# X: nObs x nEqs matrix
# return: a nObs x nEqs matrix G
mr.evalG_R <- function(y, X, beta) {
  tX <- t(X)
  # yXb <- y - c(beta %*% tX)
  yXb <- y - c(X %*% beta)
  G <- sweep(tX, MARGIN = 2, yXb, `*`)
  return(t(G))
}

# adjusted EM by chen-et-al2008
mr.evalGadj_R <- function(y,X,beta) {
  G <- mr.evalG_R(y,X,beta)
  n <- length(y)
  gbar <- 1/n*colSums(G)
  an <- max(1,0.5*log(n))
  gadd <- -an*gbar
  return(rbind(G,gadd))
}

# location-scale model
mrls.evalG_R <- function(y, X, Z, beta, gamma, sig2) {
  nObs <- nrow(X)
  nBeta <- length(beta)
  nGamma <- length(gamma)
  G <- matrix(NaN, nObs, nBeta + nGamma + 1)
  eZg <- c(exp(-Z %*% gamma)) # e^{-z'gamma}
  yXbeZg <- c((y - X %*% beta)*eZg) # (y-x'beta)e^{-z'gamma}
  yXbeZg2 <- yXbeZg * yXbeZg # (y-x'beta)^2*e^{-2z'gamma}
  G[,1:nBeta] <- yXbeZg * eZg * X
  G[,nBeta+1:nGamma] <- yXbeZg2 * Z
  # G[,nBeta+nGamma+1] <- yXbeZg2 - 1
  G[,nBeta+nGamma+1] <- 1/sig2 * yXbeZg2 - 1;
  return(G)
}

# ---- quantile regression ----

# check function
rho_alpha <- function(u, alpha) {
  u * (alpha - (u <= 0))
}

# 1st derivative of check function
phi_alpha <- function(u, alpha) {
  (u <= 0) - alpha
}

# location model
# qr.evalG is the same for non-censoring and censoring for mean regression
# X: nObs x nEqs matrix
# return: a nObs x nEqs matrix G
qr.evalG_R <- function(y, X, alpha, beta) {
  tX <- t(X) # tX is nEqs x nObs
  yXb <- y - c(beta %*% tX)
  pyXb <- phi_alpha(yXb,alpha)
  G <- sweep(tX, MARGIN = 2, pyXb, `*`)
  return(t(G))
}

# location-scale model
# TODO: only works for 1 quantile now
qrls.evalG_R <- function(y, X, Z, alpha, beta, gamma, sig2, nu) {
  nObs <- nrow(X)
  nBeta <- length(beta)
  nGamma <- length(gamma)
  G <- matrix(NaN, nObs, nBeta + nGamma + 2)
  eZg <- c(exp(-Z %*% gamma)) # e^{-z'gamma}
  yXbeZg <- c((y - X %*% beta)*eZg) # (y-x'beta)e^{-z'gamma}
  yXbeZg2 <- yXbeZg * yXbeZg # (y-x'beta)^2*e^{-2z'gamma}
  G[,1:nBeta] <- yXbeZg * eZg * X
  G[,nBeta+1:nGamma] <- yXbeZg2 * Z
  G[,nBeta+nGamma+1] <- 1/sig2 * yXbeZg2 - 1;
  G[,nBeta+nGamma+2] <- phi_alpha(yXbeZg/sqrt(sig2)-nu, alpha)
  return(G)
}
# qrls.evalG_R <- function(y, X, Z, alpha, beta, gamma, sig2, nu) {
#   nObs <- nrow(X)
#   nBeta <- length(beta)
#   nGamma <- length(gamma)
#   G <- matrix(NaN, nObs, nBeta + nGamma + 2)
#   eZg <- c(exp(-Z %*% gamma))
#   yXbeZg <- c((y - X %*% beta)*eZg/sqrt(sig2)-nu)
#   pyXbeZg <- phi_alpha(yXbeZg, alpha)
#   G[,1:nBeta] <- pyXbeZg * eZg * X # times each col of X
#   G[,nBeta+1:nGamma] <- pyXbeZg * yXbeZg * Z # times each col of Z
#   G[,nBeta+nGamma+1] <- pyXbeZg
#   G[,nBeta+nGamma+2] <- yXbeZg*yXbeZg-1
#   return(G)
# }


qr.logel_R <- function(y, X, alpha, beta, max_iter = 100, rel_tol = 1e-7) {
  G <- qr.evalG(y, X, alpha, beta)
  omegahat <- omega.hat_R(G)
  return(sum(log(omegahat)))
}

# TODO: bug
qr.post_R <- function(y, X, alpha, nsamples, nburn, betaInit, sigs) {
  betaOld <- betaInit
  betaNew <- betaOld
  # betaProp <- betaNew
  betalen <- length(betaInit)
  beta_chain <- matrix(NaN,betalen,nsamples)
  logELOld <- qr.logel_R(y, X, alpha, betaOld)
  lambdaOld <- rep(0,ncol(X))
  for (ii in (-nburn+1):nsamples) {
    for (jj in 1:betalen) {
      # satisfy <- FALSE
      betaProp <- betaOld
      # betaProp[jj] <- betaOld[jj] + sigs[jj]*rnorm(1)
      betaProp[jj] <- betaOld[jj]
      G <- qr.evalG_R(y, X, alpha, betaProp)
      lambdaOut <- lambdaNR_R(G, lambdaOld = lambdaOld)
      if (!lambdaOut$convergence) break
      # satisfy <- TRUE
      lambdaNew <- lambdaOut$lambda
      lambdaOld <- lambdaNew
      Glambda <- G %*% lambdaNew
      logomegahat <- log(1/(1-Glambda)) - log(sum(1/(1-Glambda)))
      logELProp <- sum(logomegahat)
      ratio <- exp(logELProp-logELOld)
      # message(ratio)
      u <- runif(1)
      a <- min(1,ratio)
      if (u < a) {
        betaNew <- betaProp
        betaOld <- betaNew
        logELOld <- logELProp
      }
    }
    if (ii > 0) {
      beta_chain[,ii] <- betaNew
    }
  }
  return(beta_chain)
}

# ---- posterior samplers ----

post_R <- function(Gfun, nThe, nBet, nGam,
                   y, X, Z, nsamples, nburn, 
                   ThetaInit, Sigs, RvDoMcmc, 
                   max_iter = 100, rel_tol = 1e-07) {
  nThe <- nrow(ThetaInit)
  numThe <- ncol(ThetaInit)
  Theta_chain <- matrix(NA,nrow=nThe*numThe,nsamples)
  paccept <- matrix(0,nThe,numThe)
  if (missing(RvDoMcmc)) {
    RvDoMcmc <- matrix(rep(1,nThe*numThe),nThe,numThe)
  }
  
  ThetaOld <- ThetaInit
  ThetaNew <- ThetaInit
  ThetaCur <- ThetaInit
  
  if (nThe == nBet) {
    G <- Gfun(y, X, ThetaInit)
  }
  else {
    G <- Gfun(y, X, Z, 
              ThetaInit[1:nBet,],
              ThetaInit[(nBet+1):nThe],
              ThetaInit[nThe+1,])
  }
  
  lout <- lambdaNR_R(G, max_iter, rel_tol, verbose = FALSE)
  if (!lout$convergence) message("BetaInit not valid.")
  omegas <- omega.hat_R(G)
  logelOld <- logEL_R(omegas)
  logelCur <- logelOld
  
  satisfy <- FALSE
  u <- 0
  a <- 0
  ratio <- 0
  lambda <- NaN
  
  go_next <- FALSE
  for (ii in -nburn:nsamples) {
    if (ii %% 200 == 0) message("ii = ", ii)
    go_next <- FALSE
    for (kk in 1:numThe) {
      if (go_next) break
      for (jj in 1:nThe) {
        if (RvDoMcmc[jj,kk]) {
          ThetaCur <- ThetaOld
          ThetaCur[jj,kk] <- ThetaCur[jj,kk] + Sigs[jj,kk]*rnorm(1)
          satisfy <- FALSE
          
          if (nThe == nBet) {
            G <- Gfun(y,X,ThetaCur)
          }
          else {
            G <- Gfun(y,X,Z,
                      ThetaCur[1:nBet,],
                      ThetaCur[(nBet+1):nThe],
                      ThetaCur[nThe+1,])
          }
          
          lout <- lambdaNR_R(G, max_iter, rel_tol, verbose = FALSE)
          if (!lout$convergence) {
            go_next <- TRUE
            break
          }
          else {
            satisfy <- TRUE
            lambda <- lout$lambda
          }
          
          u <- runif(1)
          logomegahat <- 1/(1 - G %*% lambda)
          logomegahat <- log(logomegahat)-log(sum(logomegahat))
          logelCur <- sum(logomegahat)
          ratio <- exp(logelCur-logelOld)
          a <- min(1.0,ratio)
          if (u < a) {
            paccept[jj,kk] <- paccept[jj,kk]+1
            ThetaNew <- ThetaCur
            ThetaOld <- ThetaCur
            logelOld <- logelCur
          }
        }
      }
    }
    if (ii > 0) {
      Theta_chain[,ii] <- as.vector(ThetaNew)
    }
  }
  paccept <- paccept/(nburn+nsamples)
  return(list(Theta_chain = Theta_chain,
              paccept = paccept))
}

mwgStep_R <- function(Gfun, nThe, nBet, nGam, y, X, Z, deltas, alpha, omegas,
                      thetaCur, logelCur, idx, mwgsd, adjust = FALSE,
                      max_iter = 100, rel_tol = 1e-07) {
  accept <- FALSE
  thetaProp <- thetaCur
  thetaProp[idx] <- thetaProp[idx] + mwgsd*rnorm(1)
  
  if (idx == nBet+nGam+1 && thetaProp[idx] < 0) {
    out <- list(accept=accept, thetaCur=thetaCur, logelCur=logelCur)
    return(out)
  }
  
  if (nThe == nBet) {
    G <- Gfun(y,X,thetaProp)
    epsilons <- evalEpsilons_R(y,X,thetaProp)
  }
  else {
    if (nThe == nBet + nGam + 1){
      G <- Gfun(y,X,Z,thetaProp[1:nBet],
                thetaProp[(nBet+1):(nBet+nGam)],
                thetaProp[nBet+nGam+1])
    }
    else {
      G <- Gfun(y,X,Z,alpha,thetaProp[1:nBet],
                thetaProp[(nBet+1):(nBet+nGam)],
                thetaProp[nBet+nGam+1],
                thetaProp[nBet+nGam+2])
    }
    epsilons <- evalEpsilonsLS_R(y,X,Z,thetaProp[1:nBet],
                                 thetaProp[(nBet+1):(nBet+nGam)],
                                 thetaProp[nGam+1])
  }
  if (adjust) {
    G <- adjG_R(G)
    weights <- evalWeights_R(c(deltas,0), omegas, c(epsilons,0))
  }
  weights <- evalWeights_R(deltas, omegas, epsilons)
  lout <- lambdaNRC_R(G,weights,max_iter,rel_tol)
  if (!lout$convergence) {
    out <- list(accept=accept, thetaCur=thetaCur, logelCur=logelCur)
    return(out)
  }
  
  omegas <- omega.hat_R(G,deltas,epsilons,adjust,max_iter,rel_tol)
  logelProp <- logEL_R(omegas,epsilons,deltas,adjust)
  u <- runif(1)
  ratio <- exp(logelProp-logelCur)
  a <- min(1.0,ratio)
  if (u < a) {
    accept <- TRUE
    thetaCur <- thetaProp
    logelCur <- logelProp
  }
  out <- list(accept=accept, thetaCur=thetaCur, logelCur=logelCur)
  return(out)
}

postCens_R <- function(Gfun, nThe, nBet, nGam,
                       y, X, Z, deltas, alpha, nsamples, nburn, 
                       thetaInit, mwgSds, rvDoMcmc, adjust = FALSE,
                       max_iter = 100, rel_tol = 1e-07) {
  # nThe <- length(thetaInit)
  paccept <- rep(0,nThe)
  theta_chain <- matrix(NA,nThe,nsamples)
  
  thetaOld <- thetaInit
  thetaCur <- thetaInit
  thetaNew <- thetaInit
  
  if (missing(rvDoMcmc)) rvDoMcmc <- rep(1,nThe)
  
  if (nThe == nBet) {
    G <- Gfun(y,X,thetaCur)
    epsilons <- evalEpsilons_R(y,X,thetaCur)
  }
  else {
    if (nThe == nBet + nGam + 1){
      G <- Gfun(y,X,Z,thetaCur[1:nBet], 
                thetaCur[(nBet+1):(nBet+nGam)], 
                thetaCur[nBet+nGam+1])
    }
    else {
      G <- Gfun(y,X,Z,alpha,thetaCur[1:nBet], thetaCur[(nBet+1):(nBet+nGam)], 
                thetaCur[nBet+nGam+1], thetaCur[nBet+nGam+2])
    }
    epsilons <- evalEpsilonsLS_R(y,X,Z,
                                 thetaCur[1:nBet],
                                 thetaCur[(nBet+1):(nBet+nGam)],
                                 thetaCur[nBet+nGam+1])
  }
  if (adjust) G <- adjG_R(G)
  omegas <- omega.hat_R(G,deltas,epsilons,adjust,max_iter,rel_tol)
  logelCur <- logEL_R(omegas,epsilons,deltas,adjust)
  
  for (ii in (-nburn+1):nsamples) {
    # if (ii %% 10 == 0) {
      message("ii = ", ii)
    # }
    for (jj in 1:nThe) {
      if (rvDoMcmc[jj]) {
        if (nThe == nBet) {
          stepout <- mwgStep_R(Gfun,nThe,nBet,nGam,y=y,X=X,deltas=deltas,
                               omegas=omegas,
                               thetaCur=thetaCur,logelCur=logelCur,idx=jj,
                               mwgsd=mwgSds[jj],adjust=adjust,
                               max_iter=max_iter,rel_tol=rel_tol)
        }
        else {
          stepout <- mwgStep_R(Gfun,nThe,nBet,nGam,y=y,X=X,Z=Z,deltas=deltas,
                               alpha=alpha, omegas=omegas,
                               thetaCur=thetaCur,logelCur=logelCur,idx=jj,
                               mwgsd=mwgSds[jj],adjust=adjust,
                               max_iter=max_iter,rel_tol=rel_tol)
        }
        paccept[jj] <- paccept[jj] + stepout$accept
        thetaCur <- stepout$thetaCur
        logelCur <- stepout$logelCur
      }
    }
    if (ii > 0) {
      theta_chain[,ii] <- thetaCur
    }
  }
  paccept <- paccept/(nsamples+nburn)
  return(list(theta_chain=theta_chain,paccept=paccept))
}

# ---- bootstrap methods with hlm ---- 
qrls_cens.boot_R <- function(y,X,Z,delta,
                             tau,beta.hat,gamma.hat,sig2.hat,nu.hat,
                             nboot=1000, max_iter=500, rel_tol=1e-6,
                             multi_cycle=FALSE) {
  p <- length(beta.hat)
  q <- length(gamma.hat)
  n <- length(y)
  beta.boot <- matrix(NA,nboot,p)
  colnames(beta.boot) <- names(beta.hat)
  gamma.boot <- matrix(NA,nboot,q)
  colnames(gamma.boot) <- names(gamma.hat)
  sig2.boot <- rep(NA,nboot)
  nu.boot <- rep(NA,nboot)
  ii <- 0
  while (ii < nboot) {
    ii <- ii+1
    if (ii %% 100 == 0) message("ii = ", ii)
    inds.boot <- sample(n,n,replace=TRUE)
    y.boot <- y[inds.boot]
    X.boot <- X[inds.boot,]
    Z.boot <- Z[inds.boot,]
    delta.boot <- deltas[inds.boot]
    hlmout <- hlm(y.boot, delta.boot, X.boot, cbind(1,Z.boot),
                  max_iter,rel_tol,multi_cycle=multi_cycle)
    if(hlmout$conv) {
      beta.boot[ii,] <- hlmout$coef$beta
      gamma.boot[ii,] <- hlmout$coef$gamma[2:(q+1)]*0.5
      sig2.boot[ii] <- exp(hlmout$coef$gamma[1])
      nu.boot[ii] <- quantile((y-X %*% beta.boot[ii,])*exp(-Z %*% gamma.boot[ii,])/sqrt(sig2.boot[ii]),tau)
    }
    else {
      ii <- ii-1
      next
    }
  }
  return(list(beta.boot=beta.boot,
              gamma.boot=gamma.boot,
              sig2.boot=sig2.boot,
              nu.boot=nu.boot))
}

# qrls_cens.boot_R <- function(y,X,Z,delta,
#                              tau,beta.hat,gamma.hat,sig2.hat,nu.hat,
#                              nboot=1000, max_iter=500, rel_tol=1e-6,
#                              multi_cycle=FALSE) {
#   p <- length(beta.hat)
#   q <- length(gamma.hat)
#   n <- length(y)
#   beta.boot <- matrix(NA,nboot,p)
#   colnames(beta.boot) <- names(beta.hat)
#   gamma.boot <- matrix(NA,nboot,q)
#   colnames(gamma.boot) <- names(gamma.hat)
#   sig2.boot <- rep(NA,nboot)
#   nu.boot <- rep(NA,nboot)
#   y.hat <- c(X %*% beta.hat + sqrt(sig2.hat)*exp(Z %*% gamma.hat))
#   eps.hat <- y - y.hat
#   ii <- 0
#   while (ii < nboot) {
#     ii <- ii+1
#     if (ii %% 20 == 0) message("ii = ", ii)
#     eps.boot <- sample(eps.hat,size=n,replace=TRUE)
#     y.boot <- y.hat + eps.boot
#     hlmout <- hlm(y.boot,delta,X,cbind(1,Z),max_iter,rel_tol,multi_cycle=multi_cycle)
#     if(hlmout$conv) {
#       beta.boot[ii,] <- hlmout$coef$beta
#       gamma.boot[ii,] <- hlmout$coef$gamma[2:(q+1)]*0.5
#       sig2.boot[ii] <- exp(hlmout$coef$gamma[1])
#       nu.boot[ii] <- quantile((y-X %*% beta.boot[ii,])*exp(-Z %*% gamma.boot[ii,])/sqrt(sig2.boot[ii]),tau)
#     }
#     else {
#       ii <- ii-1
#       next
#     }
#   }
#   return(list(beta.boot=beta.boot,
#               gamma.boot=gamma.boot,
#               sig2.boot=sig2.boot,
#               nu.boot=nu.boot))
# }

qrls_cens.bootCI_R <- function(theta.hat,theta.boot,conf=.95) {
  conf <- 1-conf
  conf <- c(conf/2, 1-conf/2)
  quantile(theta.boot - theta.hat, probs = conf) + theta.hat
  # beta.ci <- quantile(theta.hat$beta.hat - theta.boot$beta.boot, probs = conf) + theta.hat$beta.hat
  # gamma.ci <- quantile(theta.hat$gamma.hat - theta.boot$gamma.boot, probs = conf) + theta.hat$gamma.hat
  # sig2.ci <- quantile(theta.hat$sig2.hat - theta.boot$sig2.boot, probs = conf) + theta.hat$sig2.hat
  # nu.ci <- quantile(theta.hat$nu.hat - theta.boot$nu.boot, probs = conf) + theta.hat$nu.hat
  # return(list(beta.ci=beta.ci,
  #             gamma.ci=gamma.ci,
  #             sig2.ci=sig2.ci,
  #             nu.ci=nu.ci))
}
