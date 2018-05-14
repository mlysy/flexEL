# ---- generate error terms with (0,1) ----

gen_eps <- function(n, dist = "norm", df = NULL, tau) {
  
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
      nu0 <- (qlnorm(alpha,mn,sn)-m)/sqrt(v)
      return(list(eps = eps, nu0 = nu0))
    }
    else return(eps)
  }
  
  else stop("dist not included yet. ")
  
}

# # normal(0,1) error
# eps <- rnorm(n)
# nu0 <- qnorm(alpha)
# 
# # t error with mean 0 var 1
# df <- 5
# v <- df/(df-2)
# eps <- rt(n, df=df)/sqrt(v)
# nu0 <- qt(alpha,df=df)/sqrt(v)
# 
# # chi-sqr with mean 0 var 1
# df <- 3
# m <- df
# v <- 2*df
# eps <- (rchisq(n, df=df)-m)/sqrt(v)
# nu0 <- (qchisq(alpha, df)-m)/sqrt(v)
# 
# # log-normal with mean 0 var 1
# # {\displaystyle [\exp(\sigma ^{2})-1]\exp(2\mu +\sigma ^{2})}
# mn <- 0
# sn <- 1
# m <- exp(mn+sn^2/2)
# v <- (exp(sn^2)-1)*exp(2*mn+sn^2)
# eps <- (rlnorm(n,mn,sn)-m)/sqrt(v)
# nu0 <- (qlnorm(alpha,mn,sn)-m)/sqrt(v)