
#---- general helper functions ----

# generate matrix of iid N(0,1)
rMNorm <- function(n, p) {
  if(missing(p)) p <- n
  matrix(rnorm(n*p), n, p)
}

# max of min of abs and rel error
max.xdiff <- function(x) {
    xdiff <- abs(diff(x))
    max(pmin(xdiff[,1], xdiff[,2]))
}

plotEL <- function(mu.seq, logel.seq, trueval, meanobs = NA, mu.name = "param") {
    plot(mu.seq, exp(logel.seq-max(logel.seq)),
         cex=0.2, xlab = mu.name, ylab = 'log EL', type = 'l')
    abline(v = trueval, col = 'red', lty=2) # true param
    abline(v = mu.seq[which.max(logel.seq)], lty=2) # mode of EL
    if (!is.na(meanobs)) {
        abline(v = meanobs, col='blue', lty=2) # mean of observation
        legend('topleft',legend=c('true param', 'logEL mode', 'observed mean'),
               lty = c(2,2,2), col = c('red','black','blue'), cex = 0.6)
    }
    else{
        legend('topleft',legend=c('true param', 'logEL mode'),
               lty = c(2,2), col = c('red','black'), cex = 0.6)
    }
    return(mu.seq[which.max(logel.seq)]) # return the mode
}

norm_pdf <- function(ly, x) {
    py <- exp(ly - max(ly))
    py/sum(py)/(x[2]-x[1])
}

MaxRelErr <- function(lambdaNew, lambdaOld) {
    relErr <- abs((lambdaNew - lambdaOld) / (lambdaNew + lambdaOld))
    return(max(relErr))
}

#---- non-censoring EL ----

log.star <- function(x, n) {
    cond <- x >= 1/n
    ans <- rep(NA, length(x))
    ans[cond] <- log(x[cond])
    ans[!cond] <- -1/2 * n^2 * x[!cond]^2 + 2*n*x[!cond] - 3/2 - log(n)
    ans
}

log.star1 <- function(x, n) {
    cond <- x >= 1/n
    ans <- rep(NA,length(x))
    ans[cond] <- 1/(x[cond])
    ans[!cond] <- -n^2*x[!cond] + 2*n
    return(ans)
}

log.star2 <- function(x, n) {
    cond <- x >= 1/n
    ans <- rep(NA,length(x))
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
        rho <- rep(NA, nObs)
        for(jj in 1:nObs) {
            rho[jj] <- log.star1(Glambda[jj],nObs)
            Q2 <- Q2 + log.star2(Glambda[jj],nObs) * (G[,jj] %*% t(G[,jj]))
        }
        Q1 <- -G %*% rho
        lambdaNew <- lambdaOld - solve(Q2,Q1) # TODO: step size?
        maxErr <- MaxRelErr(lambdaNew, lambdaOld)
        if (verbose && (nIter %% 5 == 0)){
          message("nIter = ", nIter)
          message("err = ", maxErr)
        }
        if (maxErr < rel_tol) break;
        lambdaOld <- lambdaNew
    }
    notconv <- (ii == max_iter && maxErr > rel_tol)
    if(notconv) lambdaNew <- rep(NA, nEqs)
    # c(lambdaNew) to make sure lambda is a vector
    output <- list(lambda=c(lambdaNew), convergence=!notconv)
    return(output)
}

# 
# G is nObs x nEqs matrix
omega.hat.NC_R <- function(G, max_iter = 100, rel_tol = 1e-07, verbose = FALSE) {
    lambdaOut <- lambdaNR_R(G = G, max_iter, rel_tol, verbose)
    # nIter <- lambdaout$nIter
    # maxErr <- lambdaout$maxErr
    # notconv <- nIter == max_iter && maxErr > rel_tol # 1 if not converged
    conv <- lambdaOut$convergence # 1 if converged
    if (!conv) {
        nObs <- nrow(G)
        omegahat <- rep(0,nObs)
    }
    else {
        lambdahat <- lambdaOut$lambda
        omegahat <- c(1/(1-t(lambdahat) %*% t(G)) / sum(1/(1-t(lambdahat) %*% t(G))))
    }
    # returns a vector of omegahat
    return(omegahat)
    # return(list(omegas=omegahat, convergence=conv))
}

# G is nObs x nEqs matrix
logEL_R <- function(G, deltas, epsilons, max_iter = 100, rel_tol = 1e-7) {
  # non-censored case:
  if (missing(deltas) && missing(epsilons)) {
    omegas <- omega.hat.NC_R(G, max_iter, rel_tol)
    if (sum(omegas)== 0) return(-Inf)
    else return(sum(log(omegas)))
  }
  # censored case:
  else {
    if (length(deltas) != length(epsilons)) {
        stop("deltas and epsilons have inconsistent dimensions.")
    }
    if (nrow(G) != length(deltas)) {
      stop("deltas and G have inconsistent dimensions.")
    }
    omegas <- omega.hat.EM_R(G, deltas, epsilons, max_iter, rel_tol)
    if (sum(omegas)== 0) return(-Inf)
    else {
      epsOrd <- order(epsilons) # ascending order of epsilons
      # print(epsOrd)
      n <- length(omegas)
      psos <- rep(0,n)
      for (ii in 1:n) {
        psos[ii] <- evalPsos_R(ii, epsOrd, omegas) 
      }
      # print(psos)
      return(sum(deltas*log(omegas)+(1-deltas)*log(psos)))
    }
  }
}

# TODO: changed
mrls.logel_R <- function(y, X, Z, beta, gamma) {
    # max_iter = 100, rel_tol = 1e-07, verbose = FALSE
    G <- mrls.evalG_R(y, X, Z, beta, gamma)
    # lambda <- lambdaNR(G = G, max_iter, rel_tol, verbose)
    omegahat <- omega.hat.NC_R(G = G)
    sum(log(omegahat))
}

#---- censoring EL ----

# log.sharp and related are unique to censoring
# Note: x and q must be of the same length
log.sharp <- function(x, q) {
    cond <- x >= q
    ans <- rep(NA,length(x))
    ans[cond] <- log(x[cond])
    ans[!cond] <- -1/(2*q[!cond]^2)*x[!cond]^2 + 2/q[!cond]*x[!cond] - 3/2 + log(q[!cond])
    return(ans)
}

log.sharp1 <- function(x, q) {
    cond <- x >= q
    ans <- rep(NA,length(x))
    ans[cond] <- 1/(x[cond])
    ans[!cond] <- -1/(q[!cond]^2)*x[!cond] + 2/q[!cond]
    return(ans)
}

log.sharp2 <- function(x, q) {
    cond <- x >= q
    ans <- rep(NA,length(x))
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
    rho <- rep(NA,nObs)
    Q2 <- matrix(rep(0,nEqs*nEqs), nEqs, nEqs)
    for (jj in 1:nObs) {
      rho[jj] <- log.sharp1(Glambda[jj], weights[jj]);
      Q2 <- Q2 + weights[jj]*log.sharp2(Glambda[jj], weights[jj])*(G[,jj] %*% t(G[,jj]))
    }
    Q1 <- G %*% (rho * weights)
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
  if (notconv) lambdaNew <- rep(NA, nEqs)
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

# find a solution in the null space of G s.t. >= 0 and sum to 1. 
# G is nObs x nEqs 
# library(MASS)
# Warning: this does not guarantee to find a solution
# searchSol <- function(G) {
#     n <- nrow(G)
#     p <- ncol(G)
#     con.mat1 <- cbind(G,-G)
#     con.mat <- rbind(con.mat1, colSums(con.mat1))
#     con.dir <- c(rep(">=",n),"==")
#     con.rhs <- c(rep(0,n),1)
#     sol <- matrix(nrow=150, ncol=n) 
#     for (i in 1:n) { 
#         obj.in <- rnorm(p)
#         obj.in <- cbind(obj.in, -obj.in)
#         out <- lp (objective.in = obj.in, 
#                    const.mat = con.mat, 
#                    const.dir = con.dir, 
#                    const.rhs = con.rhs)
#         sol[i,] <- con.mat1 %*% out$solution 
#     } 
#     sol <- unique(round(sol, digits=10)) 
#     omegas <- t(sol[1,])
# }

# G is nObs x nEqs matrix 
omega.hat.EM_R <- function(G, deltas, epsilons, max_iter = 100, rel_tol = 1e-7, verbose=FALSE) {
    n <- nrow(G)
    m <- ncol(G)
    err <- Inf
    lambdaOld <- rep(0,m)
    nIter <- 0
    # initialize omegas with uncensored solution 
    omegas <- omega.hat.NC_R(G, max_iter, rel_tol, verbose)
    if (sum(omegas) == 0) return(rep(0,n))
    for (ii in 1:max_iter) {
        nIter <- ii
        # E step: calculating weights
        # weights <- getweights_R(omegas, delta)
        weights <- evalWeights_R(deltas, omegas, epsilons)
        # print("weights = ")
        # print(weights)
        # M step:
        lambdaOut <- lambdaNRC_R(G, weights, max_iter, rel_tol, verbose, lambdaOld)
        # TODO: what if not converged ?? use a random weights and continue ? 
        if (!lambdaOut$convergence) {
          message("lambdaNRC did not converge in EM")
          return(rep(0,n))
        }
        # if (!lambdaOut$convergence) {
        #     message("randomly modify omegas..")
        #     omegas <- abs(omegas + rnorm(n))
        #     omegas <- omegas/sum(omegas)
        #     lambdaNew <- lambdaOld
        #     next
        # }
        lambdaNew <- lambdaOut$lambda
        qlg <- c(sum(weights) + lambdaNew %*% t(G))
        # print(qlg)
        # print("lambdaNew = ")
        # print(lambdaNew)
        # print("lambdaOld = ")
        # print(lambdaOld)
        # if (ii == 2) break
        omegas <- weights/qlg
        # print("omegas = ")
        # print(omegas)
        omegas <- omegas/sum(omegas)
        # print("omegas = ")
        # print(omegas)
        err <- MaxRelErr(lambdaNew,lambdaOld)
        if (verbose && iter %% 20 == 0) {
            message("iter = ", iter)
            message("err = ", err)
        }
        if (err < rel_tol) break
        lambdaOld <- lambdaNew
    }
    notconv <- (nIter == max_iter && err > rel_tol) # TRUE if not converged 
    if (notconv) {
        omegas = rep(0,n)
    }
    return(omegas)
}

# wrapper function for censor and non-censor omega.hat
omega.hat_R <- function(G, deltas, epsilons, max_iter = 100, rel_tol = 1e-7, verbose = FALSE) {
    if (missing(deltas) && missing(epsilons)) {
        omegas <- omega.hat.NC_R(G, max_iter, rel_tol, verbose)
    }
    else {
        omegas <- omega.hat.EM_R(G, deltas, epsilons,
                                   max_iter, rel_tol, verbose)
    }
    return(c(omegas))
}

library(MASS) # for use of Null
# wrapper of omega.hat_R for optimCheck
# if omegas is optimal, then x == rep(1,n-p) should be optimal
omega.check <- function(x, omegas, G, deltas, epsilons) {
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
    return(sum(deltas*log(xNG)+(1-deltas)*log(psos)))
  }
}

#---- mean regression ----

# Note: returns a nObs x nEqs matrix G
# GfunCensMean <- function(y, X, beta) {
#     n <- length(y)
#     p <- length(beta)
#     if (nrow(X) != p) stop("dimensions of X and beta don't match")
#
#     G <- matrix(rep(NA,2*n),2,n)
#     for (ii in 1:n) {
#         G[1,ii] <- y[ii] - t(X[,ii]) %*% beta # mean(eps) = 0
#         G[2,ii] <- G[1,ii]^2-1 # var(eps) = 1
#     }
#     return(t(G))
# }

#  (location model)
# mr.evalG is the same for non-censoring and censoring for mean regression
# X: nObs x nEqs matrix
# return: a nObs x nEqs matrix G
mr.evalG_R <- function(y, X, beta) {
    tX <- t(X)
    yXb <- y - c(beta %*% tX)
    G <- sweep(tX, MARGIN = 2, yXb, `*`)
    return(t(G))
}

#---- quantile regression ----

rho_alpha <- function(u, alpha) {
    u * (alpha - (u <= 0))
}

phi_alpha <- function(u, alpha) {
    (u <= 0) - alpha
}

#  (location model)
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

qr.logel_R <- function(y, X, alpha, beta, max_iter = 100, rel_tol = 1e-7) {
    G <- qr.evalG(y, X, alpha, beta)
    omegahat <- omega.hat_R(G)
    return(sum(log(omegahat)))
}

qr.post_R <- function(y, X, alpha, nsamples, nburn, betaInit, sigs) {
    betaOld <- betaInit
    betaNew <- betaOld
    # betaProp <- betaNew
    betalen <- length(betaInit)
    beta_chain <- matrix(NA,betalen,nsamples)
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


## #---- Location-scale model ----

mrls.evalG_R <- function(y, X, Z, beta, gamma) {
    nObs <- nrow(X)
    nBeta <- length(beta)
    nGamma <- length(gamma)
    G <- matrix(NA, nObs, nBeta + nGamma + 1)
    yXb <- y - X %*% beta
    gZe <- exp(-2 * (Z %*% gamma))
    WW <- c(yXb * gZe)
    G[,1:nBeta] <- WW * X
    WW <- c(yXb * WW)
    G[,nBeta + 1:nGamma] <-  WW * Z
    G[,nBeta + nGamma + 1] <- WW - 1
    G
}

## LSevalG <- function(y, X, theta) {
##     nObs <- ncol(X)
##     nEqs <- nrow(X)
##     G <- matrix(rep(NA,nEqs*2*nObs), nEqs*2, nObs) # nEq*2 x nobs
##     beta <- theta[1:nEqs]
##     gamma <- theta[(nEqs+1):2*nEqs]
##     yXb <- t(y) - beta %*% X
##     gXe <- exp(-2*gamma %*% X)
##     G[1:nEqs,] <- sweep(X, MARGIN = 2, yXb*gXe)
##     yXb2 <- yXb * yXb
##     G[(nEqs+1):2*nEqs,] <- sweep(X, MARGIN = 2, yXb2*gXe)
##     return(G)
## }


