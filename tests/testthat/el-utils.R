#-------------- helper functions for EL ----------------

# interfacing C++ through R is done with the Rcpp package.
# it's a bit tricky but plenty of documentation available online.
# to handle projects with multiple C++ files (and headers) Rcpp requires you
# to build a package.  i've done this millions of times so happy to help.

## require(Rcpp)
## require(RcppEigen)
## sourceCpp(file = "NewtonRaphsonEL.cpp", verbose = TRUE)

# TODO: add comments as necessary

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

# TODO: is Q2 pd?
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

# inline double InnerEL::MaxRelErr(const Ref<const VectorXd>& lambdaNew,
#                                  const Ref<const VectorXd>& lambdaOld) {
#   relErr = ((lambdaNew - lambdaOld).array() / (lambdaNew + lambdaOld).array()).abs();
#   return(relErr.maxCoeff());
# }

MaxRelErr <- function(lambdaNew, lambdaOld) {
  relErr <- abs((lambdaNew - lambdaOld) / (lambdaNew + lambdaOld))
  return(max(relErr))
}

# Question: blockouter advantage (InnerEL) ?
# R implementation of lambdaNR
# Note: G is nEqs x nObs here
lambdaNR_R <- function(G, lambda0, max_iter=100, eps=1e-7, verbose = FALSE) {
  G <- t(G)
  nObs <- ncol(G)
  nEqs <- nrow(G)
  lambdaOld <- lambda0
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
  if(ii == max_iter & maxErr > eps) lambdaNew <- rep(NA, nEqs)
  output <- list(lambda=c(lambdaNew), maxErr=maxErr, nIter=nIter)
  return(output)
}

