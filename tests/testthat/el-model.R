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
qrls.evalG_R <- function(y, X, Z, alpha, beta, gamma, nu) {
  nObs <- nrow(X)
  nBeta <- length(beta)
  nGamma <- length(gamma)
  G <- matrix(NaN, nObs, nBeta + nGamma + 1)
  eZg <- c(exp(-Z %*% gamma))
  yXbeZg <- c((y - X %*% beta)*eZg-nu)
  pyXbeZg <- phi_alpha(yXbeZg, alpha)
  G[,1:nBeta] <- pyXbeZg * eZg * X # times each col of X
  G[,nBeta+1:nGamma] <- pyXbeZg * yXbeZg * Z # times each col of Z
  G[,nBeta+nGamma+1] <- pyXbeZg
  return(G)
}


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
                   y, X, nsamples, nburn, 
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

