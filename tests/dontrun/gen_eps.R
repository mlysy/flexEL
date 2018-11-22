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

# calculate the k-th moment of nct
# Note: k should be an integer
# nct.kmom_R <- function(k,df,ncp) {
#   if (df <= k) stop("k-th moment does not exit for df <= k.")
#   if (k == 1) {
#     return(ncp*sqrt(df/2)*gamma((df-1)/2)/gamma(df/2))
#   }
#   else if (k == 2) {
#     return(df*(1+ncp^2)/(df-2)-ncp^2*df/2*(gamma((df-1)/2)/gamma(df/2))^2)
#   }
#   else {
#     val <- (df/2)^(k/2)*gamma((df-k)/2)/gamma(df/2)*exp(-ncp^2/2)
#     encp2 <- exp(ncp^2/2)
#     if (k == 3) {
#       return(val*(3*ncp*encp2 + ncp^3*encp2))
#     }
#     else if (k == 4) {
#       return(val*(3*encp2 + 6*ncp^2*encp2 + ncp^4*encp2))
#     }
#     tltmat <- matrix(0,nrow=k+1,ncol=k) # transpose of last term matrix: coef matrix for the last term in the multiplication
#     tltmat[seq(from=2,to=k*(k+1),by=(k+2))] <- 1
#     ltmat <- t(tltmat)
#     ltmat[2,1] <- 1
#     else {
#       for (ii in 3:k) {
#         ltmat[ii,ii-1] <- ltmat[ii-1,ii-2]+ii-1
#       }
#       for (ii in 4:k) {
#         if (ii %% 2 == 0) {
#           for (jj in 1:(ii-3)) {
#             if (jj == 1) ltmat[ii,jj] <- ltmat[ii-1,jj+1]*jj
#             else ltmat[ii,jj] <- ltmat[ii-1,jj-1] + ltmat[ii-1,jj+1]*jj
#           }
#         }
#         else {
#           for (jj in 2:(ii-3)) {
#             ltmat[ii,jj] <- ltmat[ii-1,jj-1]+ltmat[ii-1,jj+1]*jj
#           }
#         }
#       }
#     }
# 
#     # TODO: return using ltmat
#   }
# }

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