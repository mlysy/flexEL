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


