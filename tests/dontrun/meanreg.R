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
