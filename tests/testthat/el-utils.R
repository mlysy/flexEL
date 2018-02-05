
#---- General helper functions ---- 

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
# Note: returns a nObs x nEqs matrix G
mr.evalG_R <- function(y, X, beta) {
    yXb <- y - beta %*% X
    G <- sweep(X, MARGIN = 2, yXb, `*`)
    return(t(G))
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

Qfun <- function(lambda, G) {
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
lambdaNR_R <- function(G, max_iter=100, eps=1e-7, verbose = FALSE) {
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
        if (maxErr < eps) break;
        lambdaOld <- lambdaNew
    }
    if(ii == max_iter && maxErr > eps) lambdaNew <- rep(NA, nEqs)
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

# function to calculate qs from ws
get_qs <- function(ws) {
    wsps <- rep(NA,n) # ws partial sum
    for (ii in 1:n) wsps[ii] <- sum(ws[ii:n])
    # W is an upper triangular matrix corresponding to w.tilde in eq (2.10)
    W <- matrix(rep(0,n*n),n,n) 
    for (jj in 1:n) W[jj,(jj:n)] <- ws[jj:n]/wsps[jj]
    qs <- rep(NA,n)
    for (kk in 1:n) {
        qs[kk] <- delta[kk] + (1-delta) %*% W[,kk]
    } 
    return(qs)
}

# for mle.check
# Note: requires G, qs to be in the environment, and G is xObs x nEqs
QfunCens <- function(lambda) {
    G <- t(G)
    qs_sum <- sum(qs)
    G_list <- split(G, rep(1:ncol(G), each = nrow(G)))
    ls <- mapply(function(gg,qq) log.sharp(x = qs_sum + sum(lambda*gg), qq),
                 G_list, qs)
    return(qs %*% ls)
}

# R implementation of lambdaNRC
# Note: input G is nObs x nEqs here
lambdaNRC_R <- function(G, lambda0, qs, maxIter = 100, eps = 1e-7, verbose = FALSE) {
    G <- t(G)
    nObs <- ncol(G) 
    nEqs <- length(lambda0) 
    lambdaOld <- lambda0
    lambdaNew <- lambdaOld
    nIter <- 0
    # newton-raphson loop
    for (ii in 1:maxIter) {
        nIter <- ii
        # Q1 and Q2
        Glambda <- lambdaOld %*% G
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
        if (maxErr < eps) {
            break;
        }
        lambdaOld <- lambdaNew # complete cycle
        if (verbose && nIter %% 20 == 0){
            message("nIter = ", nIter) 
            message("err = ", maxErr) 
        } 
    }
    output<- list(lambda = lambdaNew, nIter = nIter, maxErr = maxErr)
    return(output)
}

