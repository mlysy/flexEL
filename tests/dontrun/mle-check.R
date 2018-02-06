#' Numerical check for MLE
#'
#' Given a log-likelihood function and a potential MLE, checks whether each one-dimensional version of the log-likelihood which varies one parameter at a time with all others at the MLE is indeed maximized at the MLE value.
#' @param loglik loglikelihood function.  Takes a single vector argument.
#' @param theta.mle The potential MLE.
#' @param itheta indices of one dimensional functions to evaluate and plot.  Defaults to all parameters.
#' @param theta.names Optional vector of parameter names for plotting.
#' @param theta.rng Optional two-row matrix giving the plotting limits for each parameter.  Defaults to theta.mle +/- .5 * abs(theta.mle)
#' @param refit If \code{TRUE}, recalculates the range so that drop is more or less the same on either side of \code{theta.mle}.
#' @param layout Optional vector giving the number of rows and columns in the plot.  For \code{p} parameters, defaults to \code{c(nr, nc)}, where \code{nr = floor(p)} and \code{nc = ceiling(p/nr)}.
#' @return Invisibly returns \code{NULL}.  The purpose of this function is to plot the one-dimensional log-likelihoods.
mle.check <- function(loglik, theta.mle, itheta, theta.names, theta.rng,
                      refit = TRUE, layout, debug = FALSE) {
  ntheta <- length(theta.mle) # number of parameters
  if(missing(itheta)) itheta <- 1:ntheta
  if(is.logical(itheta)) itheta <- which(itheta) # convert T/F's to indices
  if(missing(theta.names)) {
    theta.names <- paste0("theta[",1:ntheta,"]")
    # converts to expression so symbol "theta_i" is plotted
    theta.names <- parse(text = theta.names)
  }
  if(missing(theta.rng)) {
    theta.rng <- rbind(theta.mle - .5 * abs(theta.mle),
                       theta.mle + .5 * abs(theta.mle))
  }
  # shorten theta.names and theta.rng if necessary
  ntheta2 <- length(itheta)
  if(length(theta.names) > ntheta2) theta.names <- theta.names[itheta]
  if(ncol(theta.rng) > ntheta2) theta.rng <- theta.rng[,itheta,drop=FALSE]
  # set up plot
  opar <- par(no.readonly = TRUE) # save specs of current plot
  # plot size
  if(missing(layout)) {
    layout <- floor(sqrt(ntheta2))
    layout <- c(layout, ceiling(ntheta2/layout))
  }
  # for loop for plotting
  par(mfrow = layout, mar = c(2,2.5,2.5,0), oma = c(3, 3, .5, .5))
  for(ii in 1:ntheta2) {
    ith <- itheta[ii]
    theta.seq <- seq(from = theta.rng[1,ii],
                     to = theta.rng[2,ii], len = 200)
    for(jj in 1:2) {
      # evaluate likelihood fixing all components except one
      theta.ll <- sapply(theta.seq, function(thetai) {
        theta <- theta.mle
        theta[ith] <- thetai
        loglik(theta)
      })
      if(jj == 1 && refit) {
        if(debug) browser()
        vth <- !is.na(theta.ll) & theta.ll > -Inf # valid values
        lth <- theta.seq < theta.mle[ith] # on the left of mle
        rth <- theta.seq > theta.mle[ith] # on the right
        # larger of the min value on each size
        lbd <- max(min(theta.ll[vth & lth]), min(theta.ll[vth & rth]))
        # rescale theta.seq to be on this range
        ibd <- c(which.min(ifelse(vth & lth, abs(theta.ll-lbd), Inf)),
                 which.min(ifelse(vth & rth, abs(theta.ll-lbd), Inf)))
        theta.seq <- seq(theta.seq[ibd[1]], theta.seq[ibd[2]], len = 200)
      } else break
    }
    # plot loglikelihood and add estimated value
    plot(theta.seq, theta.ll, type = "l",
         xlab = "", ylab = "")
    title(main = theta.names[ii], cex.main = 2)
    abline(v = theta.mle[ith], col = "red")
  }
  # labels in margin
  mtext(side = 2, text = "Log-Likelihood",
        line = 1, outer = TRUE)
  mtext(side = 1, text = "Parameter",
        line = 1, outer = TRUE)
  par(opar) # restore plot parameters
}
