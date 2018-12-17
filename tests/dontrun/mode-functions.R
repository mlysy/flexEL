## Note: the gradient in HMC paper only works for theta in the support of EL
## so if the initial point is far away this may not work.
## In fact, log EL may not be convex in the parameters so the gradient is not 
## guaranteed to work

# ---- shared functions ----

#'
#'@param lambda is length-\cote{p} vector
#'@param omegas is length-\code{nObs} vector
#'@param gradmat is \code{nObs x p} matrix returned by gradMat
#'@return value of the gradient of the logEL evaluated at theta
logELgrad_R <- function(omegas, lambda, gradlist) {
  grad <- sapply(gradlist, function(x) lambda %*% x)
  return(c(-grad %*% omegas))
}

#'
#'@param lambda is length-\cote{p} vector
#'@param omegas is length-\code{nObs} vector
#'@param gradmat is \code{nObs x p} matrix returned by gradMat
#'@return value of the gradient of the logEL evaluated at theta
# logELCensgrad_R <- function(omegas, deltas, epsilons, lambda, gradlist, weights) {
#   epsOrd <- order(epsilons) # ascending order of epsilons
#   n <- length(omegas)
#   psos <- rep(0,n)
#   for (ii in 1:n) {
#     psos[ii] <- evalPsos_R(ii, epsOrd, omegas) 
#   }
#   # numerical stability: watch out for extremely small negative values
#   omegas[abs(omegas) < 1e-10/length(omegas)] <- 1e-10
#   # omegas_ <- deltas*log(omegas)+(1-deltas)*log(psos)
#   grad <- sapply(gradlist, function(x) lambda %*% x)
#   # return(c(-grad %*% omegas_))
#   return(c(-grad %*% (weights*omegas)))
# }

# ---- mr functions ----

#'
#'@param X is \code{nObs x p} matrix
#'@param y is length-\code{nObs} vector
#'@param theta is length-\cote{p} vector
#'@details ...
#'@return a list of gradient values of length \code{nObs}
mr.deltaG_R <- function(y, X, theta) {
  lx <- split(X, row(X))
  dg <- lapply(lx, function(x) tcrossprod(x,x))
  return(dg)
}

#'
#'@param X is \code{nObs x p} matrix
#'@param y is length-\code{nObs} vector
#'@param beta is length-\cote{p} vector
#'@details ...
#'@return
mr.neglogEL_R <- function(y, X, beta) {
  G <- mr.evalG_R(y, X, beta)
  lout <- lambdaNR_R(G)
  if (lout$conv) {
    lambda <- lout$lambda
    omegas <- c(1/(1-t(lambda) %*% t(G)) / sum(1/(1-t(lambda) %*% t(G))))
    res <- -sum(log(omegas)) # negative logEL
    gradlist <- mr.deltaG_R(y, X, beta)
    grad <- -logELgrad_R(omegas, lambda, gradlist) # negative gradient
    attr(res, "gradient") <- grad
  }
  else {
    # TODO: if not converged, what should be the gradient...??
    res <- Inf
    attr(res, "gradient") <- rep(Inf)
  }
  return(res)
}

# ---- qr functions ----

# smoothed indicator function 1(x <= 0)
ind.smooth_R <- function(x, s=10) {
  return(1/(1+exp(s*x)))
}

# 1st derivative of ind.smooth
ind1.smooth_R <- function(x, s=10) {
  return(-s*exp(s*x)/(1+exp(s*x))^2)
}

# 2nd derivative of ind.smooth
ind2.smooth_R <- function(x, s=10) {
  return((-s*s*exp(s*x)*(1+exp(s*x))+2*s*s*exp(2*s*x))/(1+exp(s*x))^3)
}

# smoothed check function
rho.smooth_R <- function(x, tau, s=10) {
  return(x*(tau-ind.smooth_R(x, s)))
}

# 1st derivate of rho_smooth
rho1.smooth_R <- function(x, tau, s=10) {
  return(tau-ind.smooth_R(x,s)-x*ind1.smooth_R(x,s))
}

# 2nd derivative of rho_smooth
rho2.smooth_R <- function(x, tau, s=10) {
  retval <- -2 * ind1.smooth_R(x,s) - x*ind2.smooth_R(x,s)
  return(retval)
}

#'
#'@param X is \code{nObs x p} matrix
#'@param y is length-\code{nObs} vector
#'@param tau is a scalar of quantile percentage in (0,1)
#'@param beta is length-\cote{p} vector
#'@details ...
#'@return a list of gradient values of length \code{nObs}
qr.deltaG.smooth_R <- function(y, X, tau, beta, s=10) {
  lx <- split(X, row(X))
  dg <- mapply(function(x,y) {
    eps <- y - x %*% beta
    tcrossprod(rho2.smooth_R(eps,tau,s)[1]*x,x) # take [1] is because rho2.smooth_R returns a 1x1 matrix
  }, lx, y, SIMPLIFY = FALSE)
  return(dg)
}

# smoothed version of qr.evalG
qr.evalG.smooth_R <- function(y,X,tau,beta,s=10) {
  tX <- t(X)
  yXb <- y - c(beta %*% tX)
  pyXb <- rho1.smooth_R(yXb,tau,s)
  G <- sweep(tX, MARGIN = 2, pyXb, `*`)
  return(t(G))
}

#'
#'@param X is \code{nObs x p} matrix
#'@param y is length-\code{nObs} vector
#'@param beta is length-\cote{p} vector
#'@details ...
#'@return
qr.neglogEL.smooth_R <- function(y, X, tau, beta, s=10) {
  G <- qr.evalG.smooth_R(y, X, tau, beta, s)
  # lout <- lambdaNR_R(G)
  lambda <- lambdaNR(G)
  # if (lout$conv) {
  if (!anyNA(lambda)) {
    # lambda <- lout$lambda
    omegas <- c(1/(1-t(lambda) %*% t(G)) / sum(1/(1-t(lambda) %*% t(G))))
    res <- -sum(log(omegas)) # negative logEL
    gradlist <- qr.deltaG.smooth_R(y, X, tau, beta, s)
    grad <- -logELgrad_R(omegas, lambda, gradlist) # negative gradient
    attr(res, "gradient") <- grad
  }
  else {
    # TODO: if not converged, what should be the gradient...??
    message("qr.neglogEL.smooth_R: lambdaNR not coverged.")
    res <- Inf
    attr(res, "gradient") <- rep(Inf,length(beta))
  }
  return(res)
}

# ---- mrls functions ----


# ---- qrls.smooth functions ----
qrls.deltaG.smooth_R <- function(y, X, Z, tau, beta, gamma, sig2, nu, s=10) {
  lx <- split(X, row(X))
  lz <- split(Z, row(Z))
  p <- length(beta)
  q <- length(gamma)
  d <- p+q+2 # Note: this is for when the number of quantile levels is 1
  dg <- mapply(function(x,z,y) {
    # some consts to be used multiple times
    emzg <- c(exp(-z %*% gamma)) # e to the minuz z * gamma (should be a scalar)
    em2zg <- emzg * emzg # e to the minuz 2 * z * gamma
    ymxb <- c(y - x %*% beta) # y minus x * beta (should be a scalar)
    ymxb2 <- ymxb * ymxb # (y minus x * beta)^2
    sig <- sqrt(sig2) # TODO: what if negative?
    eps <- 1/sig*emzg*ymxb
    rsemm <- rho2.smooth_R(eps-nu,tau,s) # rho second derivative of eps minus nu
    
    # pre-allocate space for hessian
    mat <- matrix(NA,nrow=d,ncol=d)
    
    # derivative w.r.t beta
    mat[1:p,1:p] <- -em2zg*tcrossprod(x,x)
    mat[1:p,(p+1):(p+q)] <- 2*em2zg*ymxb*tcrossprod(x,z)
    mat[1:p,p+q+1] <- -2/sig2*em2zg*ymxb*x
    mat[1:p,d] <- -1/sig*rsemm*emzg*x
    
    # symmetric 
    mat[(p+1):(p+q),1:p] <- t(mat[1:p,(p+1):(p+q)])
    mat[p+q+1,1:p] <- t(mat[1:p,p+q+1])
    mat[d,1:p] <- t(mat[1:p,d])
    
    # derivative w.r.t. gamma
    mat[(p+1):(p+q),(p+1):(p+q)] <- 2*em2zg*ymxb2*tcrossprod(z,z)
    mat[(p+1):(p+q),p+q+1] <- -2/sig2*em2zg*ymxb2*z
    mat[(p+1):(p+q),d] <- -1/sig*rsemm*emzg*ymxb*z
    
    # symmetric
    mat[p+q+1,(p+1):(p+q)] <- t(mat[(p+1):(p+q),p+q+1])
    mat[d,(p+1):(p+q)] <- t(mat[(p+1):(p+q),d])
    
    # derivative w.r.t. sig
    mat[p+q+1,p+q+1] <- -2/(sig2*sig)*em2zg*ymxb2
    mat[p+q+1,d] <- -1/sig2*rsemm*emzg*ymxb
    
    # symmetric
    mat[d,p+q+1] <- mat[p+q+1,d]
    
    # derivative w.r.t nu
    mat[d,d] <- -rsemm
    
    # return
    mat
  }, lx, lz, y, SIMPLIFY = FALSE)
  return(dg)
}

qrls.evalG.smooth_R <- function(y, X, Z, tau, beta, gamma, sig2, nu, s=10) {
  nObs <- nrow(X)
  nBeta <- length(beta)
  nGamma <- length(gamma)
  G <- matrix(NaN, nObs, nBeta + nGamma + 2)
  eZg <- c(exp(-Z %*% gamma)) # e^{-z'gamma}
  yXbeZg <- c((y - X %*% beta)*eZg) # (y-x'beta)e^{-z'gamma}
  yXbeZg2 <- yXbeZg * yXbeZg # (y-x'beta)^2*e^{-2z'gamma}
  G[,1:nBeta] <- yXbeZg * eZg * X
  # G[,nBeta+1:nGamma] <- yXbeZg2 * Z
  # G[,nBeta+1:nGamma] <- (1-yXbeZg2) * Z
  G[,nBeta+1:nGamma] <- (1-yXbeZg2/sig2) * Z
  G[,nBeta+nGamma+1] <- 1/sig2 * yXbeZg2 - 1;
  G[,nBeta+nGamma+2] <- rho1.smooth_R(yXbeZg/sqrt(sig2)-nu, tau, s)
  return(G)
}

qrls.neglogEL.smooth_R <- function(y, X, Z, tau, theta, s=10, ret_grad=TRUE) {
  p <- ncol(X)
  q <- ncol(Z)
  beta <- theta[1:p]
  gamma <- theta[(p+1):(p+q)]
  sig2 <- theta[p+q+1]
  nu <- theta[p+q+2]
  G <- qrls.evalG.smooth_R(y, X, Z, tau, beta, gamma, sig2, nu, s)
  lambda <- lambdaNR(G)
  if (!anyNA(lambda)) {
    omegas <- omega.hat(G)
    res <- -sum(log(omegas)) # negative logEL
    if (ret_grad) {
      gradlist <- qrls.deltaG.smooth_R(y, X, Z, tau, beta, gamma, sig2, nu, s)
      grad <- -logELgrad_R(omegas, lambda, gradlist) # negative gradient
      attr(res, "gradient") <- grad
    }
  }
  else {
    # TODO: if not converged, what should be the gradient...??
    message("qrls.neglogEL_R: lambdaNR not coverged.")
    res <- rep(Inf)
    if (ret_grad) {
      attr(res, "gradient") <- rep(Inf,length(beta))
    }
  }
  return(res)
}

# ---- mr.cens functions ----
mr_cens.neglogEL.smooth_R <- function(y, X, deltas, beta, s=10) {
  G <- mr.evalG_R(y, X, beta)
  epsilons <- evalEpsilons_R(y,X,beta)
  oout <- omega.hat.EM.smooth_R(G,deltas,epsilons,s)
  if (oout$conv) {
    lambda <- oout$lambda
    omegas <- oout$omegas
    weights <- oout$weights
    res <- -logEL.smooth_R(omegas,epsilons,deltas)
    # gradlist <- mr.deltaG_R(y, X, beta)
    # grad <- -logELCensgrad_R(omegas, deltas, epsilons, lambda, gradlist, weights) # negative gradient
    # attr(res, "gradient") <- grad
    # attr(res, "gradient") <- grad(-logEL.smooth_R,)
  }
  else {
    # TODO: if not converged, what should be the gradient...??
    res <- Inf
    attr(res, "gradient") <- Inf
  }
  return(res)
}

# ---- mrls.cens functions ---- 

# ---- qrls.cens functions ----
qrls_cens.neglogEL.smooth_R <- function(y, X, Z, deltas, tau, theta, s=10) {
  nBet <- ncol(X)
  nGam <- ncol(Z)
  beta <- theta[1:nBet]
  gamma <- theta[(nBet+1):(nBet+nGam)]
  sig2 <- theta[nBet+nGam+1]
  if (sig2 < 0) return(Inf)
  nu <- theta[nBet+nGam+2]
  G <- qrls.evalG.smooth_R(y, X, Z, tau, beta, gamma, sig2, nu, s)
  if (anyNA(G)) return(Inf)
  epsilons <- evalEpsilonsLS_R(y,X,Z,beta,gamma,sig2)
  oout <- omega.hat.EM.smooth_R(G,deltas,epsilons,s)
  if (oout$conv) {
    lambda <- oout$lambda
    omegas <- oout$omegas
    weights <- oout$weights
    res <- -logEL.smooth_R(omegas,epsilons,deltas)
    # gradlist <- mr.deltaG_R(y, X, beta)
    # grad <- -logELCensgrad_R(omegas, deltas, epsilons, lambda, gradlist, weights) # negative gradient
    # attr(res, "gradient") <- grad
    # attr(res, "gradient") <- grad(-logEL.smooth_R,)
  }
  else {
    # TODO: if not converged, what should be the gradient...??
    res <- Inf
    attr(res, "gradient") <- Inf
  }
  return(res)
}
