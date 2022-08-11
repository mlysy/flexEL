## censoring example data
beta0 <- c(1,2)
n <- 200
X <- cbind(1, rnorm(n))
eps <- rnorm(n)
y <- c(X %*% beta0 + eps)


library(flexEL)
mr_neglogel <- function(beta, y, X, gel) {
  G <- flexEL::mr_evalG(y, X, beta)
  return(-gel$logel(G))
}

gel <- flexEL::GenEL$new(n_obs = n, n_eqs = 2) # initalize an GenEL object
gel$supp_adj <- TRUE # turn on support correction

nlm(mr_neglogel, c(0.75, 1.25), y, X, gel)



mr_neglogcel <- function(beta, y_obs, X, cel) {
  G <- flexEL::mr_evalG(y, X, beta)
  -gel$logel(G)
}

cel <- flexEL::CensEL$new(n_obs = n, n_eqs = 2)
cel$supp_adj <- TRUE
cel$smooth <- TRUE # turn on continuity correction
cel$smooth_s <- 1 # set tuning parameter for continuity correction

z <- rnorm(2*mean(y), n = n) # censoring variable
delta <- y <= z # censoring indicators
y_obs <- y
y_obs[!delta] <- z[!delta]

nlm(mr_neglogcel, c(0.75, 1.25), y_obs, X, cel)


# new coding 
mr_evalG_R <- function(y, X, beta) {
  tX <- t(X)
  yXb <- y - c(X %*% beta)
  G <- sweep(tX, MARGIN = 2, yXb, `*`)
  return(t(G))
}


mr_neglogEL_R <- function(y, X, beta) {
  G <- mr_evalG_R(y, X, beta) # G matrix based on your regression problem
  return(-flexEL::logEL(G = G, supp_adj = FALSE)) # support correction not on
}


gen_nct_eps <- function(n, df, ncp) {
  m <- ncp*sqrt(df/2)*gamma((df-1)/2)/gamma(df/2)
  v <- df*(1+ncp^2)/(df-2)-ncp^2*df/2*(gamma((df-1)/2)/gamma(df/2))^2
  eps <- (rt(n, df = df, ncp = ncp)-m)/sqrt(v)
}
n <- 500 # number of observations
b <- c(0.5,1) # beta_0 = 0.5, beta_1 = 1
eps <- gen_nct_eps(n, df = 20, ncp = 1) # a re-centered right-skewed non-central t distribution
X <- cbind(1, rnorm(n)) # n x 2 covariate matrix (intercept included)
y <- X %*% b + eps # length n response vector



beta_init <- coef(lm(y ~ X-1)) # obtain initial value
nlmout <- nlm(f = mr_neglogEL_R, p = beta_init, y = y, X = X)
nlmout$estimate


mr_dGdb_R <- function(y, X, beta) {
  lx <- split(X, row(X))
  dg <- lapply(lx, function(x) tcrossprod(x,x))
  return(dg)
}



mr_dGdb_R(y,X,beta_init)

beta <- beta_init

mr_neglogEL_R(y,X,beta)
beta <- c(1.4491182,1.9705781)


mr_neglogEL_R <- function(y, X, beta) {
  G <- mr_evalG_R(y, X, beta)
  res_lst <- flexEL::logEL(G, supp_adj = FALSE, grad = TRUE) 
  neglogel <- -res_lst$logel
  dldG <- -res_lst$grad
  dGdb <- mr_dGdb_R(y, X, beta)
  grad_mat <- matrix(NA, nrow = nrow(dldG), ncol = ncol(dldG))
  for (ii in 1:ncol(dldG)) {
    grad_mat[,ii] <- dGdb[[ii]] %*% dldG[ii,]
  }
  grad <- colSums(grad_mat)
  attr(neglogel, "gradient") <- grad
  return(neglogel)
}


nlmout <- nlm(f = mr_neglogEL_R, p = beta_init, y = y, X = X)
nlmout$estimate





## the example in the jss paper

require(flexEL)
require(R6)

meanreg_logel <- function(y, X, beta) {
  yXb <- y - X %*% beta
  G <- as.numeric(yXb) * X
  flexEL::logEL(G = G, supp_adj = TRUE)
}

MeanReg <- R6Class(
  classname = "MeanReg",
  inherit = GenEL,
  
  private = list(
    .X = NULL,
    .y = NULL
  ),
  
  public = list(
    
    #' Class constructor.
    initialize = function(y, X) {
      n_obs <- nrow(X)
      n_eqs <- ncol(X)
      super$initialize(n_obs, n_eqs)
      private$.X <- X
      private$.y <- y
    },
    
    #' Calculate G matrix.
    evalG = function(beta) {
      yXb <- private$.y - private$.X %*% beta
      as.numeric(yXb) * private$.X
    },
    
    #' Calculate the log-EL.
    logel = function(beta) {
      G <- self$evalG(beta)
      super$logel(G)
    }
  )
)


beta0 <- c(1,2)
n <- 200
X <- cbind(1, rnorm(n))
eps <- rnorm(n)
y <- c(X %*% beta0 + eps)



mr_neglogel <- function(beta, y, X, gel) {
  G <- flexEL::mr_evalG(y, X, beta)
  return(-gel$logel(G))
}

gel <- flexEL::GenEL$new(n_obs = n, n_eqs = 2) # initalize an GenEL object
gel$supp_adj <- TRUE # turn on support correction

nlm(mr_neglogel, c(0.75, 1.25), y, X, gel)

mrel <- MeanReg$new(y, X)
mrel$supp_adj <- TRUE
mrel$logel(c(0.75, 1.25))


meanreg_logel(y,X,c(0.75, 1.25))






MeanReg <- R6Class(
  classname = "MeanReg",
  inherit = GenEL,
  
  private = list(
    .X = NULL,
    .y = NULL
  ),
  
  public = list(
    
    #' Class constructor.
    initialize = function(y, X) {
      n_obs <- nrow(X)
      n_eqs <- ncol(X)
      super$initialize(n_obs, n_eqs)
      private$.X <- X
      private$.y <- y
    },
    
    #' Calculate G matrix.
    evalG = function(beta) {
      yXb <- private$.y - private$.X %*% beta
      as.numeric(yXb) * private$.X
    },
    
    #' Calculate the log-EL.
    logel = function(beta) {
      G <- self$evalG(beta)
      super$logel(G)
    },
    
    #' Calculate the gradient of logel w.r.t. G 
    dldG = function(beta) {
      G <- self$evalG(beta)
     super$logel_grad(G)$grad
    },
    
    
    #' Calculate the gradient of logel w.r.t. beta
    dldb = function(beta){
      # construct dGdb
      lx <- split(private$.X, row(private$.X))
      # calculate dGdb
      dGdb <- lapply(lx, function(x) tcrossprod(x,x))
      dldG_mat <- self$dldG(beta)
      grad_mat <- matrix(NA, nrow = nrow(dldG_mat), ncol = ncol(dldG_mat))
      for (ii in 1:ncol(dldG_mat)) {
        grad_mat[,ii] <- dGdb[[ii]] %*% dldG_mat[ii,]
      }
      grad <- colSums(grad_mat)
      grad
      
    }
    
  )
)


mrel <- MeanReg$new(y, X)
mrel$supp_adj <- TRUE
mrel$logel(c(0.75, 1.25))
mrel$dldb(c(1, 2))

dGdb <- mr_dGdb_R(y, X, beta)

beta <- c(1,2)
G <- mr_evalG_R(y, X, beta)
res_lst <- flexEL::logEL(G, supp_adj = TRUE, grad = TRUE) 
neglogel <- -res_lst$logel
dldG <- -res_lst$grad
dGdb <- mr_dGdb_R(y, X, beta)
grad_mat <- matrix(NA, nrow = nrow(dldG), ncol = ncol(dldG))
for (ii in 1:ncol(dldG)) {
  grad_mat[,ii] <- dGdb[[ii]] %*% dldG[ii,]
}
grad <- colSums(grad_mat)




mr_dGdb_R <- function(y, X, beta) {
  lx <- split(X, row(X))
  dg <- lapply(lx, function(x) tcrossprod(x,x))
  return(dg)
}

mr_neglogEL_R <- function(y, X, beta) {
  G <- mr_evalG_R(y, X, beta)
  res_lst <- flexEL::logEL(G, supp_adj = TRUE, grad = TRUE) 
  neglogel <- -res_lst$logel
  dldG <- -res_lst$grad
  dGdb <- mr_dGdb_R(y, X, beta)
  grad_mat <- matrix(NA, nrow = nrow(dldG), ncol = ncol(dldG))
  for (ii in 1:ncol(dldG)) {
    grad_mat[,ii] <- dGdb[[ii]] %*% dldG[ii,]
  }
  grad <- colSums(grad_mat)
  attr(neglogel, "gradient") <- grad
  return(neglogel)
}


mr_neglogEL_R(y,X,c(1, 2))




## in the jss paper
beta0 <- c(1,2)
n <- 200
X <- cbind(1, rnorm(n))
eps <- rnorm(n)
y <- c(X %*% beta0 + eps)


require(flexEL)
#' create G matrix 
mr_evalG_R <- function(y, X, beta) {
  yXb <- y - X %*% beta
  G <- as.numeric(yXb) * X
  return(G)
}

#' return logEL value with the gel class
mr_neglogel <- function(beta, y, X, gel) {
  G <- flexEL::mr_evalG(y, X, beta)
  return(-gel$logel(G))
}

gel <- flexEL::GenEL$new(n_obs = n, n_eqs = 2) # initialize an GenEL object
gel$supp_adj <- TRUE # turn on support correction

nlm(mr_neglogel, c(0.75, 1.25), y, X, gel)


require(R6)

meanreg_logel <- function(y, X, beta) {
  yXb <- y - X %*% beta
  G <- as.numeric(yXb) * X
  flexEL::logEL(G = G, supp_adj = TRUE)
}

MeanReg <- R6Class(
  classname = "MeanReg",
  inherit = GenEL,
  
  private = list(
    .X = NULL,
    .y = NULL
  ),
  
  public = list(
    
    #' Class constructor.
    initialize = function(y, X) {
      n_obs <- nrow(X)
      n_eqs <- ncol(X)
      super$initialize(n_obs, n_eqs)
      private$.X <- X
      private$.y <- y
    },
    
    #' Calculate G matrix.
    evalG = function(beta) {
      yXb <- private$.y - private$.X %*% beta
      as.numeric(yXb) * private$.X
    },
    
    #' Calculate the log-EL.
    logel = function(beta) {
      G <- self$evalG(beta)
      super$logel(G)
    },
    
    #' Calculate the gradient of logel w.r.t. G 
    dldG = function(beta) {
      G <- self$evalG(beta)
      super$logel_grad(G)$grad
    },
    
    
    #' Calculate the gradient of logel w.r.t. beta
    dldb = function(beta){
      # construct dGdb
      lx <- split(private$.X, row(private$.X))
      # calculate dGdb
      dGdb <- lapply(lx, function(x) tcrossprod(x,x))
      # calculate dldG
      dldG_mat <- self$dldG(beta)
      grad_mat <- matrix(NA, nrow = nrow(dldG_mat), ncol = ncol(dldG_mat))
      for (ii in 1:ncol(dldG_mat)) {
        grad_mat[,ii] <- dGdb[[ii]] %*% dldG_mat[ii,]
      }
      grad <- colSums(grad_mat)
      grad
    }
  )
)



#' check the logEL is the same 
#' functional result
meanreg_logel(y, X, c(0.75, 1.25))

#' object-oriented result
mrel <- MeanReg$new(y, X)
mrel$supp_adj <- TRUE
mrel$logel(c(0.75, 1.25))


#' check the gradient is the same 
#' functional result 
#' the function returning dGdb
mr_dGdb_R <- function(y, X, beta) {
  lx <- split(X, row(X))
  dg <- lapply(lx, function(x) tcrossprod(x,x))
  return(dg)
}

#' the function returning dldb
mr_neglogEL_R <- function(y, X, beta) {
  G <- mr_evalG_R(y, X, beta)
  res_lst <- flexEL::logEL(G, supp_adj = TRUE, grad = TRUE) 
  neglogel <- -res_lst$logel
  dldG <- -res_lst$grad
  dGdb <- mr_dGdb_R(y, X, beta)
  grad_mat <- matrix(NA, nrow = nrow(dldG), ncol = ncol(dldG))
  for (ii in 1:ncol(dldG)) {
    grad_mat[,ii] <- dGdb[[ii]] %*% dldG[ii,]
  }
  grad <- colSums(grad_mat)
  attr(neglogel, "gradient") <- grad
  return(neglogel)
}
mr_neglogEL_R(y,X,c(1, 2))
attr(mr_neglogEL_R(y,X,c(1, 2)),"gradient")
#' object-oriented result 
mrel$dldb(c(1, 2))


identical(attr(mr_neglogEL_R(y,X,c(1, 2)),"gradient"),-mrel$dldb(c(1, 2)))



# quantile example 
n <- 20
p <- 2
X <- replicate(p-1, rnorm(n))
X <- cbind(1, X)
beta0 <- c(1, 2)
y <- c(X %*% beta0) + rnorm(n)



alpha <- 0.5
beta <- c(1, 2)
qr_evalG(y, X, alpha, beta, s = 1)


## simulate some data ##
n <- 20
p <- 2
q <- 2
X <- replicate(p-1, rnorm(n))
X <- cbind(1, X) # X may have an intercept term, if so, it should be explicit
Z <- replicate(q, rnorm(n)) # Z shall not have an intercept term
beta0 <- c(1, 2)
gamma0 <- c(0.5, 0.25)
sig20 <- 0.5
y <- c(X %*% beta0 + sqrt(sig20)*exp(Z %*% gamma0)*rnorm(n)) # with N(0,1) error term

## calculate G matrix given data and certain parameter values ##
# a single quantile level (with continuity correction)
alpha <- 0.5
beta <- c(1, 2)
gamma <- c(0.5, 0.25)
sig2 <- 0.5
sigma <- sqrt(sig2)
nu <- 0.2
G <- qrls_evalG(y, X, Z, alpha, beta, gamma, sig2, nu, s = 1)


yXb <- as.numeric(y - X %*% beta)
eZg <- exp(as.numeric(Z %*% gamma))

# beta
G[,1:2]
(yXb/eZg^2) * X   
# gamma     
G[,3:4]      
(1 - (yXb/(sigma*eZg))^2) * Z
# sigma
G[,5]
(yXb/(sigma*eZg))^2 - 1     
# nu
check_smooth_dx <- function(x, tau, s) (tau - 1/(1+exp(s*x))) + x*(s*exp(s*x)/((1+exp(s*x))^2))
  
G[,6]
check_smooth_dx(yXb/(sigma*eZg) - nu, tau = .5, s = 1)


#' Quantile regression G function.
quantreg_evalG <- function(y, X, Z, beta, gamma, sigma, nu, tau, s) {
  yXb <- as.numeric(y - X %*% beta)
  eZg <- exp(as.numeric(Z %*% gamma))
  G <- cbind((yXb/eZg^2) * X,                 # beta
             (1 - (yXb/(sigma*eZg))^2) * Z,   # gamma
             (yXb/(sigma*eZg))^2 - 1,         # sigma
             check_smooth_dx(yXb/(sigma*eZg) - nu, tau, s))
  G
}


require(R6)
QuantReg <- R6Class(
  classname = "QuantReg",
  inherit = CensEL,
  
  private = list(
    .X = NULL,
    .Z = NULL,
    .u = NULL,
    .delta = NULL,
    .s = NULL,
    .tau = NULL
  ),
  
  public = list(
    
    initialize = function(u, delta, X, Z) {
      n_obs <- nrow(X)
      n_eqs <- ncol(X) + ncol(Z) + 2
      super$initialize(n_obs, n_eqs)
      private$.X <- X
      private$.Z <- Z
      private$.u <- u
      private$.delta <- delta
      private$.s <- s
      private$.tau <- tau
    },
    
    evalG = function(beta, gamma, sigma, nu) {
      quantreg_evalG(y = private$.u,
                     X = private$.X,
                     Z = private$.Z,
                     beta = beta,
                     gamma = gamma,
                     sigma = sigma,
                     nu = nu,
                     tau = private$.tau,
                     s = private$.s)
    },
    
    resid = function(beta, gamma, sigma) {
      uXb <- private$.u - as.numeric(private$.X %*% beta)
      eZg <- sigma * exp(as.numeric(private$.Z %*% gamma))
      uXb/eZg
    },
    
    logel = function(beta, gamma, sigma, nu) {
      eps <- self$resid(beta, gamma, sigma)
      G <- evalG(beta, gamma, sigma, nu, private$.tau, private$.s)
      super$logel(G, private$.delta, eps)
    }
  )
)

# CensEL$new(n_obs, n_eqs)
flexEL::logEL(G, supp_adj = FALSE, grad = TRUE) 

qrls_neglogel <- function(theta){
   beta <- theta[1:2] 
   gamma <- theta[3:4]
   sig2 <- theta[5]
   nu <- theta[6]
   alpha <- .5
   G <- qrls_evalG(y, X, Z, alpha, beta, gamma, sig2, nu, s = 1)
  -flexEL::logEL(G, supp_adj = FALSE)
   
}
theta <- c(1,2.1,.9,.25,.4,.5)
qrls_neglogel(theta)

optim(theta,qrls_neglogel, method = "BFGS")


# multiple quantile levels (with continuity correction)
alpha <- c(0.25, 0.75)
beta <- c(1, 2)
gamma <- c(0.5, 0.25)
sig2 <- 0.5
nu <- c(-0.5, 0.5)
qrls_evalG(y, X, Z, alpha, beta, gamma, sig2, nu, s = 1)




# check the gradient of check_smooth_dx


#' Smoothed check function.
#' #' Smoothed indicator function.
indicator_smooth <- function(x, s) 1/(1 + exp(s * x))

check_smooth <- function(x, tau, s) {
  x * (tau - indicator_smooth(x, s))
  }

#' Derivative of smooth checked function wrt `x`
check_smooth_dx <- function(x, tau, s) {
  (tau - 1/(1+exp(s*x))) + x*(s*exp(s*x)/((1+exp(s*x))^2))
}


check_smooth_dx(x = 2, tau = .5, s = 1)

check_smooth <- function(x, tau, sm) {
  x * (tau - 1/(1 + exp(sm * x)))
}
numDeriv::grad(check_smooth,  x = 1, sm = 1, tau = .5)

numDeriv::grad(sin, pi)

x <- rnorm(20)
numDeriv::grad(check_smooth,  x,tau = 1,s = .5)
check_smooth_dx(x = x, tau = .5, s = 1)

numDeriv::grad(check_smooth, x) - check_smooth_dx(x = x, tau = .5, s = 1)



func0 <- function(x,y,s){ sin(x) +y + x*s}
numDeriv::grad(func0, x =1, y = 9,  s = 2)



func0 <- function(x,s,tau = .5){ 
  x * (tau - 1/(1 + exp(s * x)))
}
numDeriv::grad(func0, x =1,s = .5)



check_smooth <- function(x, tau, sm) {
  x * (tau - 1/(1 + exp(sm * x)))
}
numDeriv::grad(check_smooth,  x = 1, sm = 1, tau = .5)



# the quantile regression example



## simulate some data ##
require(flexEL)
# n is the sample size
n <- 100
# p is the length of beta
p <- 2
# q is the length of gamma 
q <- 1
# X is the location covariate variable
X <- replicate(p-1, rnorm(n))
X <- cbind(1, X) # X may have an intercept term, if so, it should be explicit
# Z is the scale covariate variable 
Z <- replicate(q, rnorm(n)) # Z shall not have an intercept term
# true value for beta
beta0 <- c(.5, 1)
# true value for gamma
gamma0 <- -.5
# true value for sigma
sig20 <- 1
# the error term eps ~ N(0,1)
eps <- rnorm(n)
# simulate the y data with error term eps ~ N(0,1)
y <- c(X %*% beta0 + sqrt(sig20)*exp(Z %*% gamma0)*eps) # with N(0,1) error term

## calculate G matrix given data and certain parameter values ##
# a single quantile level (with continuity correction) with 0.75 quantile
alpha <- 0.75 
# initial value 
beta <- c(.5, 1) 
gamma <- -.5 
sig2 <- 1 
sigma <- sqrt(sig2)
nu <- .6 


## censored data

c <- rnorm(n = n, mean = 1.35, sd = 1) # censoring variable
cap_eps <- c <= eps # censoring indicators
y_allcen <- c(X %*% beta0 + sqrt(sig20)*exp(Z %*% gamma0)*c)
y_obs <- y
y_obs[cap_eps] <- y_allcen[cap_eps]
delta <- as.numeric(cap_eps)
delta_rep <- -1*(delta -1 )


G <- flexEL::qrls_evalG(y, X, Z, alpha, beta, gamma, sig2, nu, s = 1)


# CensEL$new(n_obs, n_eqs)
flexEL::logEL(G, supp_adj = FALSE) 

qrls_neglogel <- function(theta,alpha){
  beta <- theta[1:2] 
  gamma <- theta[3]
  sig2 <- theta[4]
  nu <- theta[5]
  G <- qrls_evalG(y, X, Z, alpha, beta, gamma, sig2, nu, s = 1)
  -flexEL::logEL(G, supp_adj = FALSE)
  
}

optim(theta_init,qrls_neglogel,alpha = .75)$par


#### censoring data
# object-oriented result
mrel <- MeanReg$new(y, X)
mrel$supp_adj <- TRUE
# the logEL value in the OOP implementation
mrel$logel(c(0.75, 1.25))



cel <- flexEL::CensEL$new(n_obs = n, n_eqs = 2)
cel$supp_adj <- TRUE
cel$smooth <- TRUE # turn on continuity correction
cel$smooth_s <- 1 # set tuning parameter for continuity correction

c <- rnorm(n = n, mean = 1.35, sd = 1) # censoring variable
cap_eps <- c <= eps # censoring indicators
y_allcen <- c(X %*% beta0 + sqrt(sig20)*exp(Z %*% gamma0)*c)
y_obs <- y
y_obs[cap_eps] <- y_allcen[cap_eps]

### check the censoring is correct 
# 1 cap_eps  = FALSE
cap_eps[1]
# it means  c > eps
# so it should use eps 
X[1,] %*% beta0 + sqrt(sig20)*exp(Z[1] %*% gamma0)*eps[1]
X[1,] %*% beta0 + sqrt(sig20)*exp(Z[1] %*% gamma0)*c[1]
y_obs[1]
# 29 cap_eps = TRUE
cap_eps[29]
# it means  c > eps
# so it should use eps 
X[29,] %*% beta0 + sqrt(sig20)*exp(Z[29] %*% gamma0)*eps[29]
X[29,] %*% beta0 + sqrt(sig20)*exp(Z[29] %*% gamma0)*c[29]
y_obs[29]


# initialize the quanreg class
#--- Helper functions ----------------------------------------------------------

#' Smoothed indicator function.
indicator_smooth <- function(x, s) 1/(1 + exp(s * x))

#' Smoothed check function.
check_smooth <- function(x, tau, s) x * (tau - indicator_smooth(x, s))

#' Derivative of smooth checked function wrt `x`.
check_smooth_dx <- function(x, tau, s) (tau - 1/(1+exp(s*x))) + x*(s*exp(s*x)/((1+exp(s*x))^2))


#' Quantile regression G function.
quantreg_evalG <- function(y, X, Z, beta, gamma, sigma, nu, tau, s) {
  yXb <- as.numeric(y - X %*% beta)
  eZg <- exp(as.numeric(Z %*% gamma))
  G <- cbind((yXb/eZg^2) * X,                 # beta
             (1 - (yXb/(sigma*eZg))^2) * Z,   # gamma
             (yXb/(sigma*eZg))^2 - 1,         # sigma
             check_smooth_dx(yXb/(sigma*eZg) - nu, tau, s))
  G
}

QuantReg <- R6Class(
  classname = "QuantReg",
  inherit = CensEL,
  
  private = list(
    .X = NULL,
    .Z = NULL,
    .u = NULL,
    .delta = NULL,
    .s = NULL,
    .tau = NULL
  ),
  
  public = list(
    
    initialize = function(u, delta, X, Z, s, tau) {
      n_obs <- nrow(X)
      n_eqs <- ncol(X) + ncol(Z) + 2
      super$initialize(n_obs, n_eqs)
      private$.X <- X
      private$.Z <- Z
      private$.u <- u
      private$.delta <- delta
      private$.s <- s
      private$.tau <- tau
    },
    
    evalG = function(beta, gamma, sigma, nu) {
      quantreg_evalG(y = private$.u,
                     X = private$.X,
                     Z = private$.Z,
                     beta = beta,
                     gamma = gamma,
                     sigma = sigma,
                     nu = nu,
                     tau = private$.tau,
                     s = private$.s)
    },
    
    resid = function(beta, gamma, sigma) {
      uXb <- private$.u - as.numeric(private$.X %*% beta)
      eZg <- sigma * exp(as.numeric(private$.Z %*% gamma))
      uXb/eZg
    },
    
    logel = function(beta, gamma, sigma, nu) {
      eps <- self$resid(beta, gamma, sigma)
      G <- self$evalG(beta, gamma, sigma, nu)
      super$logel(G, private$.delta, eps)
    }
  )
)
## simulate some data ##
require(flexEL)
# n is the sample size
n <- 500
# p is the length of beta
p <- 2
# q is the length of gamma 
q <- 1
# X is the location covariate variable
X <- replicate(p-1, rnorm(n))
X <- cbind(1, X) # X may have an intercept term, if so, it should be explicit
# Z is the scale covariate variable 
Z <- replicate(q, rnorm(n)) # Z shall not have an intercept term
# true value for beta
beta0 <- c(.5, 1)
# true value for gamma
gamma0 <- -.5
# true value for sigma
sig20 <- 1
# the error term eps ~ N(0,1)
eps <- rnorm(n)
# simulate the y data with error term eps ~ N(0,1)
y <- c(X %*% beta0 + sqrt(sig20)*exp(Z %*% gamma0)*eps) # with N(0,1) error term

## calculate G matrix given data and certain parameter values ##
# a single quantile level (with continuity correction) with 0.75 quantile
alpha <- 0.75 
# initial value 
beta <- c(.5, 1) 
gamma <- -.5 
sig2 <- 1 
sigma <- sqrt(sig2)
nu <- .6 


## censored data

c <- rnorm(n = n, mean = 1.35, sd = 1) # censoring variable
cap_eps <- c <= eps # censoring indicators
y_allcen <- c(X %*% beta0 + sqrt(sig20)*exp(Z %*% gamma0)*c)
y_obs <- y
y_obs[cap_eps] <- y_allcen[cap_eps]
delta <- as.numeric(cap_eps)
delta_rep <- -1*(delta -1 )

s = 1
tau = .75
delta <- as.numeric(cap_eps)
delta_rep <- -1*(delta -1 )
qrel <- QuantReg$new(y_obs, delta_rep, X, Z, s, tau)
qrel$supp_adj <- FALSE
qrel$smooth <- TRUE
qrel$smooth_s <- 1
qrel$evalG(beta = c(1,2),
           gamma = -.5,
           sigma = 1,
           nu = .6)

qrel$resid(beta = c(1,2),
           gamma = -.5,
           sigma = 1)

qrel$logel(beta = c(.5,1),
           gamma = -.5,
           sigma = 1,
           nu = .7)


theta_init
neg_quantreg_logel <- function(theta) {
  beta <- theta[1:2]
  gamma <- theta[3]
  sigma <- theta[4]
  nu <- theta[5]
  -qrel$logel(beta,gamma,sigma,nu)
}

neg_quantreg_logel(theta_init)
nlmresult <- nlm(neg_quantreg_logel, theta_init)
nlmresult

optim( theta_init,neg_quantreg_logel)
# optimout <- optim(
#   theta_init, 
#   neg_quantreg_logel,
#   method = "L-BFGS-B",
#   lower = c(-Inf,-Inf,-Inf,0,0),
#   upper = c(Inf,Inf,Inf, Inf, 1)
# )



n_obs <- nrow(X)
n_eqs <- ncol(X) + ncol(Z) + 2


G <- quantreg_evalG(y = y_obs,
                    X = X,
                    Z = Z,
                    beta = beta,
                    gamma = gamma,
                    sigma = sigma,
                    nu = nu,
                    tau = tau,
                    s = s)

uXb <- y_obs - as.numeric(X %*% beta)
eZg <- sigma * exp(as.numeric(Z %*% gamma))



epsilon <- uXb/eZg
cel <- CensEL$new(n_obs, n_eqs)
cel$smooth <- TRUE # turn on continuity correction
cel$smooth_s <- 1 # set tuning parameter for continuity correction
cel$supp_adj <- TRUE
CensEL$logel(G, delta, epsilon)


s = 1
tau = .75
delta <- as.numeric(cap_eps)
qrel <- QuantReg$new(y_obs, delta, X, Z, s, tau)
qrel$supp_adj <- FALSE
qrel$smooth <- TRUE
G_class <- qrel$evalG(beta = c(.5,1),
           gamma = -.5,
           sigma = 1,
           nu = .6)
# check the G's, which are the same now 
G_class <- qrel$evalG(beta = c(.5,1),
                      gamma = -.5,
                      sigma = 1,
                      nu = .6)
G_func <- quantreg_evalG(y = y_obs,
                    X = X,
                    Z = Z,
                    beta = beta,
                    gamma = gamma,
                    sigma = sigma,
                    nu = nu,
                    tau = tau,
                    s = s)
G_class - G_func

# check the residuals, which are the same now
uXb <- y_obs - as.numeric(X %*% beta)
eZg <- sigma * exp(as.numeric(Z %*% gamma))



epsilon_func <- uXb/eZg 
epsilon_class <- qrel$resid(beta = c(.5,1),
           gamma = -.5,
           sigma = 1)

epsilon_func - epsilon_class

# check the logEL
cel <- CensEL$new(n_obs, n_eqs)
cel$smooth <- TRUE # turn on continuity correction
cel$smooth_s <- 1 # set tuning parameter for continuity correction
cel$supp_adj <- FALSE

delta_rep <- -1*(delta -1 )
cel$logel(G_func, delta_rep, epsilon_func)

cel <- CensEL$new(n_obs, n_eqs)
cel$smooth <- TRUE # turn on continuity correction
cel$smooth_s <- 1 # set tuning parameter for continuity correction
cel$supp_adj <- FALSE


logel_reg <- function(theta){
  beta <- theta[1:2]
  gamma <- theta[3]
  sigma <- theta[4]
  nu <- theta[5]
  G_func <- quantreg_evalG(y = y_obs,
                           X = X,
                           Z = Z,
                           beta = beta,
                           gamma = gamma,
                           sigma = sigma,
                           nu = nu,
                           tau = tau,
                           s = s)
  uXb <- y_obs - as.numeric(X %*% beta)
  eZg <- sigma * exp(as.numeric(Z %*% gamma))
  epsilon_func <- uXb/eZg 
  cel <- CensEL$new(n_obs, n_eqs)
  cel$smooth <- TRUE # turn on continuity correction
  cel$smooth_s <- 1 # set tuning parameter for continuity correction
  cel$supp_adj <- FALSE
  -cel$logel(G_func, delta_rep, epsilon_func)
  
}
theta <- theta_init
logel_reg(theta_init)
logel_reg(theta_init)
optim(theta_init, logel_reg)











## wrap-up coding for jss paper


#--- Helper functions ----------------------------------------------------------

#' Smoothed indicator function.
indicator_smooth <- function(x, s) 1/(1 + exp(s * x))

#' Smoothed check function.
check_smooth <- function(x, tau, s) x * (tau - indicator_smooth(x, s))

#' Derivative of smooth checked function wrt `x`.
check_smooth_dx <- function(x, tau, s) (tau - 1/(1+exp(s*x))) + x*(s*exp(s*x)/((1+exp(s*x))^2))


#' Quantile regression G function.
quantreg_evalG <- function(y, X, Z, beta, gamma, sigma, nu, tau, s) {
  yXb <- as.numeric(y - X %*% beta)
  eZg <- exp(as.numeric(Z %*% gamma))
  G <- cbind((yXb/eZg^2) * X,                 # beta
             (1 - (yXb/(sigma*eZg))^2) * Z,   # gamma
             (yXb/(sigma*eZg))^2 - 1,         # sigma
             check_smooth_dx(yXb/(sigma*eZg) - nu, tau, s))
  G
}

QuantReg <- R6Class(
  classname = "QuantReg",
  inherit = CensEL,
  
  private = list(
    .X = NULL,
    .Z = NULL,
    .u = NULL,
    .delta = NULL,
    .s = NULL,
    .tau = NULL
  ),
  
  public = list(
    
    initialize = function(u, delta, X, Z, s, tau) {
      n_obs <- nrow(X)
      n_eqs <- ncol(X) + ncol(Z) + 2
      super$initialize(n_obs, n_eqs)
      private$.X <- X
      private$.Z <- Z
      private$.u <- u
      private$.delta <- delta
      private$.s <- s
      private$.tau <- tau
    },
    
    evalG = function(beta, gamma, sigma, nu) {
      quantreg_evalG(y = private$.u,
                     X = private$.X,
                     Z = private$.Z,
                     beta = beta,
                     gamma = gamma,
                     sigma = sigma,
                     nu = nu,
                     tau = private$.tau,
                     s = private$.s)
    },
    
    resid = function(beta, gamma, sigma) {
      uXb <- private$.u - as.numeric(private$.X %*% beta)
      eZg <- sigma * exp(as.numeric(private$.Z %*% gamma))
      uXb/eZg
    },
    
    logel = function(beta, gamma, sigma, nu) {
      eps <- self$resid(beta, gamma, sigma)
      G <- self$evalG(beta, gamma, sigma, nu)
      super$logel(G, private$.delta, eps)
    }
  )
)
## simulate some data ##
require(flexEL)
# n is the sample size
n <- 200
# p is the length of beta
p <- 2
# q is the length of gamma 
q <- 1
# X is the location covariate variable
X <- replicate(p-1, rnorm(n))
X <- cbind(1, X) # X may have an intercept term, if so, it should be explicit
# Z is the scale covariate variable 
Z <- replicate(q, rnorm(n)) # Z shall not have an intercept term
# true value for beta
beta0 <- c(.5, 1)
# true value for gamma
gamma0 <- -.5
# true value for sigma
sig20 <- 1
# the error term eps ~ N(0,1)
eps <- rnorm(n)
# simulate the y data with error term eps ~ N(0,1)
y <- c(X %*% beta0 + sqrt(sig20)*exp(Z %*% gamma0)*eps) # with N(0,1) error term

## calculate G matrix given data and certain parameter values ##
# a single quantile level (with continuity correction) with 0.75 quantile
alpha <- 0.75 
# initial value 
beta <- c(.5, 1) 
gamma <- -.5 
sig2 <- 1 
sigma <- sqrt(sig2)
nu <- qnorm(.75)
theta_init <- c(beta, gamma, sigma, nu)

## censored data

c <- rnorm(n = n, mean = 1.35, sd = 1) # censoring variable
cap_eps <- c >= eps # censoring indicators
y_allcen <- c(X %*% beta0 + sqrt(sig20)*exp(Z %*% gamma0)*c)
y_obs <- y
y_obs[cap_eps] <- y_allcen[cap_eps]
delta <- as.numeric(cap_eps)

s = 1
tau = .75
qrel <- QuantReg$new(y_obs, delta, X, Z, s, tau)
qrel$supp_adj <- FALSE
qrel$smooth <- TRUE
qrel$smooth_s <- 1
qrel$evalG(beta = c(1,2),
           gamma = -.5,
           sigma = 1,
           nu = .6)

qrel$resid(beta = c(1,2),
           gamma = -.5,
           sigma = 1)

qrel$logel(beta = c(.5,1),
           gamma = -.5,
           sigma = 1,
           nu = .7)



neg_quantreg_logel <- function(theta) {
  beta <- theta[1:2]
  gamma <- theta[3]
  sigma <- theta[4]
  nu <- theta[5]
  -qrel$logel(beta,gamma,sigma,nu)
}
theta_init
neg_quantreg_logel(theta_init)

optim( theta_init,neg_quantreg_logel)
