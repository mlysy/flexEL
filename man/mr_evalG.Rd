% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mr_evalG.R
\name{mr_evalG}
\alias{mr_evalG}
\title{Evaluate the G matrix for a mean regression model.}
\usage{
mr_evalG(y, X, beta)
}
\arguments{
\item{y}{A numeric vector of resoponses of length \code{n_obs} where \code{n_obs}
is the number of observations.}

\item{X}{A numeric matrix of covariates of dimension \code{n_obs} x \code{n_bet}
where \code{n_obs} is the number of observations and \code{b_bet} is the number
of coefficients (length of \code{beta}).}

\item{beta}{A numeric vector of coefficients.}
}
\value{
A numeric matrix of dimension \code{n_obs} x \code{n_bet}.
}
\description{
Evaluate the G matrix for a mean regression model.
}
\details{
Assuming data were generated from\preformatted{y_i = x_i'beta + eps_i, for i = 1, ..., n,
}

where \code{eps_i}'s are ~iid \code{F(eps)}, with \code{E[eps] = 0} and \code{Var[eps] = 1}.
The estimating equations for \code{beta} can be derived based on the least squares method,
that is,\preformatted{min_beta (y - X*beta)'(y - X*beta)
}

The above minimization problem is equivalent to setting it's first derivative
w.r.t \code{beta} to 0, i.e.,\preformatted{X'*(y - X*beta) = 0
}

The left-hand-side of the above equation is the \code{G} matrix returned by this function.
}
\examples{
## simulate some data ##
n <- 20
p <- 2
X <- replicate(p-1, rnorm(n))
X <- cbind(1, X)
beta0 <- c(1, 2)
y <- c(X \%*\% beta0) + rnorm(n) # with N(0,1) error term

## calculate G matrix given data and certain parameter values ##
beta <- c(1, 2)
mr_evalG(y, X, beta)
}
