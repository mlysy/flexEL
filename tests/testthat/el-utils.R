
#---- General helper functions ----

# generate matrix of iid N(0,1)
rMNorm <- function(n, p) {
  if(missing(p)) p <- n
  matrix(rnorm(n*p), n, p)
}

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
# Note: returns a nObs x nEqs matrix G
mr.evalG_R <- function(y, X, beta) {
    tX <- t(X)
    yXb <- y - c(beta %*% tX)
    G <- sweep(tX, MARGIN = 2, yXb, `*`)
    return(t(G))
}

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

mrls.logel_R <- function(y, X, Z, beta, gamma,
                         max_iter = 100, rel_tol = 1e-07, verbose = FALSE) {
  G <- mrls.evalG_R(y, X, Z, beta, gamma)
  lambda <- lambdaNR(G = G, max_iter, rel_tol, verbose)
  sum(logomegahat(lambda, t(G)))
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

# Is Q2 pd?
# solveV <- function(V, x, ldV = FALSE) {
#   C <- chol(V) # cholesky decomposition
#   if(missing(x)) x <- diag(nrow(V))
#   # solve is O(ncol(C) * ncol(x)) with triangular matrices
#   # using backward subsitution
#   ans <- backsolve(r = C, x = backsolve(r = C, x = x, transpose = TRUE))
#   if(ldV) {
#     ldV <- 2 * sum(log(diag(C)))
#     ans <- list(y = ans, ldV = ldV)
#   }
#   ans
# }

#---- Non-censoring case ----

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

# Note: input G is nObs x nEqs
Qfun <- function(lambda, G) {
    G <- t(G)
    N <- ncol(G) # nObs
    sum(apply(G, 2, function(gg) log.star(x = 1 - sum(lambda*gg), n = N)))
}

logomegahat <- function(lambdahat, G) {
    # returns a vector of log omegahat
    log(1/(1-t(lambdahat) %*% G) / sum(1/(1-t(lambdahat) %*% G)))
}

rho_alpha <- function(u, alpha) {
    u * (alpha - (u <= 0))
}

phi_alpha <- function(u, alpha) {
    (u <= 0) - alpha
}

# R implementation of lambdaNRC
# Note: input G is nObs x nEqs here
lambdaNR_R <- function(G, max_iter=100, rel_tol=1e-7, verbose = FALSE) {
    G <- t(G)
    nObs <- ncol(G)
    nEqs <- nrow(G)
    lambdaOld <- rep(0,nEqs)
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
        lambdaNew <- lambdaOld - solve(Q2,Q1)
        maxErr <- MaxRelErr(lambdaNew, lambdaOld)
        if (maxErr < rel_tol) break;
        lambdaOld <- lambdaNew
    }
    if(ii == max_iter && maxErr > rel_tol) lambdaNew <- rep(NA, nEqs)
    # c(lambdaNew) to make sure lambda is a vector
    output <- list(lambda=c(lambdaNew), maxErr=maxErr, nIter=nIter)
    return(output)
}

#---- Censoring case ----

# log.sharp and related are unique to censoring
log.sharp <- function(x, q) {
    cond <- x >= q
    ans <- rep(NA,length(x))
    ans[cond] <- log(x[cond])
    ans[!cond] <- -1/(2*q[!cond]^2)*x[!cond]^2 + 2/q[!cond]*x[!cond]
    - 3/2 + log(q[!cond])
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
# Note: requires qs to be in the environment, and G is xObs x nEqs
QfunCens <- function(lambda, G) {
    G <- t(G)
    qs_sum <- sum(qs)
    G_list <- split(G, rep(1:ncol(G), each = nrow(G)))
    ls <- mapply(function(gg,qq) log.sharp(x = qs_sum + sum(lambda*gg), qq),
                 G_list, qs)
    return(qs %*% ls)
}

# R implementation of lambdaNRC
# Note: input G is nObs x nEqs here
lambdaNRC_R <- function(G, qs, maxIter = 100, rel_tol = 1e-7, verbose = FALSE,
                        lambdaOld = NULL) {
    G <- t(G)
    nObs <- ncol(G)
    nEqs <- nrow(G)
    if (is.null(lambdaOld)) lambdaOld <- rep(0,nEqs)
    lambdaNew <- lambdaOld
    nIter <- 0
    # newton-raphson loop
    for (ii in 1:maxIter) {
        nIter <- ii
        # Q1 and Q2
        Glambda <- t(lambdaOld) %*% G
        Glambda <- sum(qs) + Glambda
        Q1 <- rep(0,nEqs)
        Q2 <- matrix(rep(0,nEqs*nEqs), nEqs, nEqs)
        for (jj in 1:nObs) {
            Q1 <- Q1 + qs[jj]*log.sharp1(Glambda[jj],qs[jj])*G[,jj]
            Q2 <- Q2 + qs[jj]*log.sharp2(Glambda[jj],qs[jj])*(G[,jj] %*% t(G[,jj]))
        }
        # update lambda
        # tryCatch(solve(Q2,Q1), finally = print(lambdaOld))
        # Problem lambdaOld flies off to +Inf, but because Q2 tends to 0
        lambdaNew <- lambdaOld - solve(Q2,Q1)
        maxErr <- MaxRelErr(lambdaNew, lambdaOld) # maximum relative error
        if (maxErr < rel_tol) {
            break;
        }
        lambdaOld <- lambdaNew # complete cycle
        if (verbose && nIter %% 20 == 0){
            message("nIter = ", nIter)
            message("err = ", maxErr)
        }
    }
    output<- list(lambda = c(lambdaNew), nIter = nIter, maxErr = maxErr)
    return(output)
}


# getqs_R <- function(ws, delta) {
#     n <- length(ws)
#     wsps <- rep(NA,n) # ws partial sum
#     for (ii in 1:n) wsps[ii] <- sum(ws[ii:n])
#     # W is an upper triangular matrix corresponding to w.tilde in eq (2.10)
#     W <- matrix(rep(0,n*n),n,n)
#     for (jj in 1:n) W[jj,(jj:n)] <- ws[jj:n]/wsps[jj]
#     # print(W)
#     qs <- rep(NA,n)
#     for (kk in 1:n) {
#         qs[kk] <- delta[kk] + (1-delta) %*% W[,kk]
#     }
#     return(qs)
# }


evalPsos_R <- function(ii, epsOrd, omegas) {
    nObs <- length(omegas)
    psos <- 0
    for (jj in 1:nObs) {
        kk <- epsOrd[jj]
        psos <- psos + omegas[kk] # sum over the omegas of eps <= than ii-th eps
        if (kk == ii) break
    }
    return(psos)
}

evalWeights_R <- function(y, X, deltas, omegas, beta) {
    nObs <- length(y)
    epsilons <- y - c(X %*% beta)
    epsOrd <- order(epsilons, decreasing = TRUE)
    psots <- rep(0,nObs)
    for (ii in 1:nObs) {
        for (jj in nObs:1) {
            kk <- epsOrd[jj]
            if (deltas[kk] == 0) {
                psots[ii] <- psots[ii] + omegas[ii]/evalPsos_R(kk, epsOrd, omegas)
            }
            if (kk == ii) break
        }
    }
    weights <- deltas + psots
    return(weights)
}

# EMEL_R <- function(G, delta, ws0, max_iter = 100, rel_tol = 1e-7, verbose=FALSE) {
EMEL_R <- function(y, X, G, delta, ws0, max_iter = 100, rel_tol = 1e-7, verbose=FALSE) {
    G <- t(G)
    n <- ncol(G)
    m <- nrow(G)
    ws <- ws0
    err <- Inf
    lambdaOld <- rep(0,m)
    nIter <- 0
    for (ii in 1:max_iter) {
        nIter <- ii
        # E step: calculating qs
        # qs <- getqs_R(ws, delta)
        qs <- evalWeights_R(y, X, deltas, omegas, beta)
        print(qs)
        # M step:
        lambdaNew <- lambdaNRC_R(t(G), qs, max_iter, rel_tol, verbose, lambdaOld)$lambda
        print(lambdaNew)
        qs_sum <- sum(qs)
        qlg <- qs_sum + lambdaNew %*% G
        ws <- qs/qlg
        ws <- ws/sum(ws)
        err <- MaxRelErr(lambdaNew,lambdaOld)
        if (verbose && iter %% 20 == 0) {
            message("iter = ", iter)
            message("err = ", err)
        }
        if (err < rel_tol) break
        lambdaOld <- lambdaNew
    }
    if (nIter == max_iter && err > rel_tol) {
        ws = rep(NA,n)
    }
    output <- list(ws = ws, lambda = lambdaNew)
    return(output)
}

## #---- Location-scale model ----
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


