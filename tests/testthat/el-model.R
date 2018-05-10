# ---- mean regression ----

# location model
# mr.evalG is the same for non-censoring and censoring for mean regression
# X: nObs x nEqs matrix
# return: a nObs x nEqs matrix G
mr.evalG_R <- function(y, X, beta) {
  tX <- t(X)
  yXb <- y - c(beta %*% tX)
  G <- sweep(tX, MARGIN = 2, yXb, `*`)
  return(t(G))
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

