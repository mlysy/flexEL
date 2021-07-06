#' @title Quantile regression EL class
#'
#' @description R6 class for EL with mean regression moment specification.
#'
#' @export
QuantRegEL <- R6::R6Class(
  
  classname = "QuantRegEL",
  
  private = list(
    .el_obj = NULL,
    .is_cens = FALSE,
    .is_ls = FALSE,
    .y = NULL,
    .X = NULL,
    .Z = NULL,
    .alpha = NULL, # quantile level(s)
    .n_obs = NULL,
    .n_bet = NULL,
    .n_gam = NULL,
    .n_eqs = NULL,
    .delta = NULL,
    .s = NULL # smooth parameter
  ),
  
  public = list(
    
    initialize = function(y, X, Z = NULL, delta = NULL, alpha, s = 0, el_opts = NULL) {
      private$.alpha = alpha
      private$.s = s
      if (length(y) != nrow(X)) {
        stop("Number of observations in y must equal to that in X.")
      }
      private$.y <- y
      private$.X <- X
      private$.n_obs <- length(y)
      private$.n_bet <- ncol(X)
      if (!is.null(Z)) {
        if (length(y) != nrow(Z)) {
          stop("Number of observations in y must equal to that in Z.")
        }
        private$.is_ls <- TRUE
        private$.Z <- Z
        private$.n_gam <- ncol(Z)
        private$.n_eqs <- private$.n_bet + private$.n_gam + 1
      }
      else {
        private$.n_eqs <- private$.n_bet
      }
      if (!is.null(delta)) {
        private$.is_cens <- TRUE
        private$.delta <- delta
      }
      if(private$.is_cens) {
        private$.el_obj <- CensEL$new(n_obs = private$.n_obs, 
                                      n_eqs = private$.n_eqs)
      } else {
        private$.el_obj <- GenEL$new(n_obs = private$.n_obs, 
                                     n_eqs = private$.n_eqs)
      }
      if (!is.null(el_opts)) {
        do.call(private$.el_obj$set_opts, el_opts)
      }
    },
    
    eval_G = function(theta) {
      
      if (length(theta) != private$.n_eqs) {
        stop("Length of theta must equal to number of estimating equations.")
      }
      
      if (!private$.is_ls) {
        qr_evalG(private$.y, private$.X, private$.alpha, theta, private$.s)
      }
      else {
        beta <- theta[1:private$.n_bet]
        gamma <- theta[(private$.n_bet+1):(private$.n_bet+private$.n_gam)]
        sig2 <- theta[private$.n_bet+private$.n_gam+1]
        if (sig2 <= 0) {
          stop("Scale parameter must be a positive number.")
        }
        qrls_evalG(private$.y, private$.X, private$.Z, 
                   beta, gamma, sig2, private$.alpha, private$.s)
      }
    },
    
    logel = function(theta, grad_out = TRUE, negate = TRUE) {
      
      # calculate G matrix 
      G <- self$eval_G(theta)
      
      # calculate logel
      if (!private$.is_cens) {
        lel <- private$.el_obj$logel(G)
      }
      else {
        if (!private$.is_ls) {
          epsilon <- c(private$.y - private$.X %*% theta)
        }
        else {
          beta <- theta[1:private$.n_bet]
          gamma <- theta[(private$.n_bet+1):(private$.n_bet+private$.n_gam)]
          sig2 <- theta[private$.n_bet+private$.n_gam+1]
          epsilon <- c((private$.y - private$.X %*% beta)/(sqrt(sig2)*exp(private$.Z %*% gamma)))
        }
        lel <- private$.el_obj$logel(G, private$.delta, epsilon)
      }
      if (negate) lel <- -lel
      
      # calculate dldt if needed
      if (!private$.is_cens && grad_out) {
        dldt <- self$logel_grad(theta, logel_out = FALSE, negate = negate)
        attr(lel, "gradient") <- dldt
      }
      else if (private$.is_cens && grad_out) {
        warning("Gradient of logel under right-censoring is not supported.")
      }
      return(lel)
    },
    
    logel_grad = function(theta, logel_out = FALSE, negate = TRUE) {
      
      # calculate G matrix 
      G <- self$eval_G(theta)
      
      # calculate dldG
      if (!private$.is_cens) {
        dldG <- private$.el_obj$logel_grad(G)$dldG
      }
      else {
        stop("Gradient of logel under right-censoring is not supported.")
      }
      
      # calculate dGdt
      if (!private$.is_ls) {
        dGdt <- qr_dGdt(private$.y, private$.X, theta, private$.alpha, private$.s)
      }
      else {
        beta <- theta[1:private$.n_bet]
        gamma <- theta[(private$.n_bet+1):(private$.n_bet+private$.n_gam)]
        sig2 <- theta[private$.n_bet+private$.n_gam+1]
        dGdt <- qrls_dGdt(private$.y, private$.X, private$.Z, 
                          beta, gamma, sig2)
      }
      
      # calculate dldt
      grad_mat <- matrix(NA, nrow = nrow(dldG), ncol = ncol(dldG))
      for (ii in seq_along(1:nrow(dldG))) {
        grad_mat[ii,] <- dGdt[[ii]] %*% dldG[ii,]
      }
      dldt <- colSums(grad_mat)
      if (negate) dldt <- -dldt
      
      # calculate logel if needed
      if (logel_out) {
        logel <- self$logel(theta, grad_out = FALSE, negate = negate)
        attr(dldt, "logel") <- logel
      }
      return(dldt)
    },
    
    get_el_obj = function() {
      return(private$.el_obj)
    }
  )
)

# n <- 100
# p <- 2
# X <- replicate(p, rnorm(n))
# X[,1] <- rep(1,n)
# beta <- rnorm(p)
# alpha <- runif(1)
# y <- c(X %*% beta) + rnorm(n) # with N(0,1) error term
# sp <- sample(1:100, 1)
# qr_el <- flexEL::QuantRegEL$new(y, X, alpha = 0.5, s = 1)
# qr_el$logel(beta)
# numDeriv::grad(qr_el$logel, beta)
# nlm(qr_el$logel, beta)
