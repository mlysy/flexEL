library(flexEL)

# ---- generate error terms with (0,1) ----

gen_eps <- function(n, dist = "norm", df = NULL, ncp = NULL, tau) {
  
  # TODO: input check here
  
  if (dist == "norm") {
    eps <- rnorm(n)
    if (!missing(tau)) {
      nu0 <- qnorm(tau)
      return(list(eps = eps, nu0 = nu0))
    }
    else return(eps)
  }
  
  else if (dist == "t") {
    v <- df/(df-2)
    eps <- rt(n, df=df)/sqrt(v)
    if (!missing(tau)) {
      nu0 <- qt(tau,df=df)/sqrt(v)
      return(list(eps = eps, nu0 = nu0))
    }
    else return(eps)
  }
  
  else if (dist == "nct") {
    if (df <= 2) stop("variance do not exist for df <= 2.")
    m <- ncp*sqrt(df/2)*gamma((df-1)/2)/gamma(df/2)
    v <- df*(1+ncp^2)/(df-2)-ncp^2*df/2*(gamma((df-1)/2)/gamma(df/2))^2
    eps <- (rt(n, df = df, ncp = ncp)-m)/sqrt(v)
    if (!missing(tau)) {
      nu0 <- (qt(tau,df=df,ncp=ncp)-m)/sqrt(v)
      return(list(eps = eps, nu0 = nu0))
    }
    else return(eps)
  }
  
  else if (dist == "chisq") {
    m <- df
    v <- 2*df
    eps <- (rchisq(n, df=df)-m)/sqrt(v)
    if (!missing(tau)) {
      nu0 <- (qchisq(tau, df)-m)/sqrt(v)
      return(list(eps = eps, nu0 = nu0))
    }
    else return(eps)
  }
  
  else if (dist == "lnorm") {
    mn <- 0
    sn <- 1
    m <- exp(mn+sn^2/2)
    v <- (exp(sn^2)-1)*exp(2*mn+sn^2)
    eps <- (rlnorm(n,mn,sn)-m)/sqrt(v)
    if (!missing(tau)) {
      nu0 <- (qlnorm(tau,mn,sn)-m)/sqrt(v)
      return(list(eps = eps, nu0 = nu0))
    }
    else return(eps)
  }
  
  # else if (dist == "lgam") {
  #   ...
  # }
  
  else stop("dist not included yet. ")
  
}

# ---- two-param mean regression ----

mr_evalG_R <- function(y, X, beta) {
  tX <- t(X)
  yXb <- y - c(X %*% beta)
  G <- sweep(tX, MARGIN = 2, yXb, `*`)
  return(t(G))
}

mr_dGdb_R <- function(y, X, beta) {
  lx <- split(X, row(X))
  dg <- lapply(lx, function(x) -tcrossprod(x,x))
  return(dg)
}

mr_neglogel_with_grad_R <- function(beta, gel, y, X) {
  G <- mr_evalG_R(y, X, beta)
  logel_lst <- gel$logel_grad(G)
  neglogel <- -logel_lst$logel
  dldG <- -logel_lst$dldG
  dGdb <- mr_dGdb_R(y, X, beta)
  grad_mat <- matrix(NA, nrow = nrow(dldG), ncol = ncol(dldG))
  for (ii in 1:nrow(dldG)) {
    grad_mat[ii,] <- dGdb[[ii]] %*% dldG[ii,]
  }
  grad <- colSums(grad_mat)
  attr(neglogel, "gradient") <- grad
  return(neglogel)
}

mr_neglogel_R <- function(beta, gel, y, X) {
  G <- mr_evalG_R(y, X, beta)
  return(-gel$logel(G))
}

mr_neglogel_grad_R <- function(beta, gel, y, X) {
  G <- mr_evalG_R(y, X, beta)
  logel_lst <- gel$logel_grad(G)
  dldG <- -logel_lst$dldG
  dGdb <- mr_dGdb_R(y, X, beta)
  grad_mat <- matrix(NA, nrow = nrow(dldG), ncol = ncol(dldG))
  for (ii in 1:nrow(dldG)) {
    grad_mat[ii,] <- dGdb[[ii]] %*% dldG[ii,]
  }
  colSums(grad_mat)
}

n <- 500
X <- cbind(1, rnorm(n))
beta0 <- c(1, 2)
eps <- gen_eps(n, dist = "chisq", df = 5, tau = 0.75)$eps
y <- c(X %*% beta0) + eps
gel <- GenEL$new(n_obs = n, n_eqs = 2)

# neglogel_mr_R(c(1, 2), gel, y, X)

b0 <- seq(from = 0, to = 2, length.out = 100)
b1 <- seq(from = 1, to = 3, length.out = 100)

ll <- rep(NA, 100)
for (ii in seq_along(b0)) {
  ll[ii] <- mr_neglogel_R(c(b0[ii], 2), gel, y, X)
}
plot(b0, ll, type = 'l')

ll <- rep(NA, 100)
for (ii in seq_along(b0)) {
  ll[ii] <- mr_neglogel_R(c(1, b1[ii]), gel, y, X)
}
plot(b1, ll, type = 'l')

beta_init <- coef(lm(y ~ X - 1))
nlm(f = mr_neglogel_with_grad_R, p = c(0.5,3), gel, y, X)

optim(c(0.5,3), fn = mr_neglogel_R, gr = mr_neglogel_grad_R, gel, y, X, method = "BFGS")
