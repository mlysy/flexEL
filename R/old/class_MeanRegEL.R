#' @title Mean regression EL class
#'
#' @description R6 class for EL with location or location-scale mean regression moment specification.
#'
#' @export
MeanRegEL <- R6::R6Class(
  
  classname = "MeanRegEL",
  
  private = list(
    .el_obj = NULL,
    .is_cens = FALSE,
    .is_ls = FALSE,
    .y = NULL,
    .X = NULL,
    .Z = NULL,
    .n_obs = NULL,
    .n_bet = NULL,
    .n_gam = NULL,
    .n_eqs = NULL,
    .delta = NULL
  ),
  
  public = list(
    
    #' @description Create a new `MeanRegEL` object. This class can be used for 
    #'   either non-censored or right-censored data, and location or location 
    #'   mean regression models.
    #' @template args-y_X_Z
    #' @template arg-delta
    #' @param el_opts A list of options for the `GenEL` or `CensEL` object. 
    #'   See `set_opts` method in the corresponding class.
    #' @return A `MeanRegEL` object.
    #' @seealso [flexEL::GenEL$set_opts()], [flexEL::CensEL$set_opts()]
    initialize = function(y, X, Z = NULL, delta = NULL, el_opts = NULL) {
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
    
    #' @description Evaluate the G matrix for a given parameter value.
    #' @param theta Length-\code{n_eqs} vector of parameters for a location or location-scale model. This vector should correspond to `beta` for a location regression model, or concatenated `beta`, `gamma`, and `sig2` for a location-scale regression model, in this order.
    #' @return A numeric matrix of dimension `n_obs x n_eqs`.
    #' @seealso [flexEL::mr_evalG()], [flexEL::mrls_evalG()]
    eval_G = function(theta) {
      
      if (length(theta) != private$.n_eqs) {
        stop("Length of theta must equal to number of estimating equations.")
      }
      
      if (!private$.is_ls) {
        mr_evalG(private$.y, private$.X, theta)
      }
      else {
        beta <- theta[1:private$.n_bet]
        gamma <- theta[(private$.n_bet+1):(private$.n_bet+private$.n_gam)]
        sig2 <- theta[private$.n_bet+private$.n_gam+1]
        if (sig2 <= 0) {
          stop("Scale parameter must be a positive number.")
        }
        mrls_evalG(private$.y, private$.X, private$.Z, 
                   beta, gamma, sig2)
      }
    },
    
    #' @description Calculate the log empirical likelihood base on the given 
    #'   parameter value and in the non-censoring case can also return the 
    #'   gradient.
    #' @param theta Length-\code{n_eqs} vector of parameters for a location or location-scale model. This vector should correspond to `beta` for a location regression model, or concatenated `beta`, `gamma`, and `sig2` for a location-scale regression model, in this order.
    #' @param grad_out If `TRUE`, add gradient as an attribute to the return value.
    #' @param negate If `TRUE`, calculate the negative log likelihood.
    #' @return A scalar with a possible attribute.
    #' @seealso [flexEL::GenEL$logel()], [flexEL::CensEL$logel()], 
    #'   [flexEL::GenEL$logel_grad()], [flexEL::CensEL$logel_grad()]
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
    
    #' @description Calculate the gradient of log empirical likelihood at the 
    #'   given parameter value for non-censoring data, can also return the 
    #'   log empirical likelihood.
    #' @param theta Length-\code{n_eqs} vector of parameters for a location or location-scale model. This vector should correspond to `beta` for a location regression model, or concatenated `beta`, `gamma`, and `sig2` for a location-scale regression model, in this order.
    #' @param logel_out If `TRUE`, add log empirical likelihood as an attribute 
    #'   to the return value.
    #' @param negate If `TRUE`, calculate the gradient of negative log likelihood.
    #' @seealso [flexEL::GenEL$logel()], [flexEL::CensEL$logel()], 
    #'   [flexEL::GenEL$logel_grad()], [flexEL::CensEL$logel_grad()]
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
        dGdt <- mr_dGdt(private$.y, private$.X, theta)
      }
      else {
        beta <- theta[1:private$.n_bet]
        gamma <- theta[(private$.n_bet+1):(private$.n_bet+private$.n_gam)]
        sig2 <- theta[private$.n_bet+private$.n_gam+1]
        dGdt <- mrls_dGdt(private$.y, private$.X, private$.Z, 
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
    
    #' @description Returns the `GenEL` or `CensEL` object for inspection or 
    #'   debugging purposes.
    #' @return A `GenEL` or `CensEL` object.
    get_el_obj = function() {
      return(private$.el_obj)
    }
  )
)