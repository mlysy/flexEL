## ---- echo=TRUE, eval=TRUE----------------------------------------------------
mr_evalG_R <- function(y, X, beta) {
  tX <- t(X)
  yXb <- y - c(X %*% beta)
  G <- sweep(tX, MARGIN = 2, yXb, `*`)
  return(t(G))
}

## ---- echo=TRUE, eval=TRUE----------------------------------------------------
mr_neglogEL_R <- function(y, X, beta) {
  G <- mr_evalG_R(y, X, beta) # G matrix based on your regression problem
  return(-flexEL::logEL(G = G, supp_adj = FALSE)) # support correction not on
}

## ---- echo=TRUE, eval=TRUE----------------------------------------------------
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

## ---- echo=TRUE, eval=TRUE----------------------------------------------------
beta_init <- coef(lm(y ~ X-1)) # obtain initial value
nlmout <- nlm(f = mr_neglogEL_R, p = beta_init, y = y, X = X)
nlmout$estimate

## ---- echo=TRUE, eval=TRUE----------------------------------------------------
mr_dGdb_R <- function(y, X, beta) {
  lx <- split(X, row(X))
  dg <- lapply(lx, function(x) tcrossprod(x,x))
  return(dg)
}

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

## ---- echo=TRUE, eval=TRUE----------------------------------------------------
nlmout <- nlm(f = mr_neglogEL_R, p = beta_init, y = y, X = X)
nlmout$estimate

## ---- echo = FALSE, results = "asis"------------------------------------------
cat("```c", readLines("example.cpp"), "```", sep = "\n")

## ---- warn=FALSE, include=TRUE, eval=TRUE-------------------------------------
Rcpp::sourceCpp("example.cpp")
bb <- c(1,2)
n_obs <- 200
n_eqs <- 2
X <- cbind(1, rnorm(n_obs))
eps <- rnorm(n_obs)
y <- X %*% bb + eps
example_logel(c(0.75, 1.25), X, y)

