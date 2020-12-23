#' R6 class for EL with general moment specification.
#' 
#' A general EL object.
GenEL <- R6::R6Class(
  
  classname = "GenEL",
  
  private = list(
    
    .n_obs = NULL,
    .n_eqs = NULL,
    # .G = NULL,
    .max_iter = 100,
    .rel_tol = 1e-7,
    .supp_corr = TRUE,
    .lambda0 = 0,
    
    #' @description Check the dimension of G is what is expected.
    #' @param G A numeric matrix.
    check_G = function(G) {
      if (nrow(G) != private$.n_obs | ncol(G) != private$.n_eqs) {
        stop("Dimension of G matrix must be `n_obs` x `n_eqs`.")
      }
    }
  ),
  
  active = list(
    
    #' @description Access or reset the value of .G.
    #' @param value Missing or a numeric matrix.
    .G = function(value) {
      if (missing(value)) private$.G
      else {
        private$check_G(value)
        private$.G <- value
      }
    }
  ),
  
  public = list(
    
    #' @description Create a new GenEL object.
    #' @template arg-G
    initialize = function(n_obs, n_eqs) {
      private$.n_obs <- n_obs
      private$.n_eqs <- n_eqs
      self$.G <- matrix(NA, nrow = n_obs, ncol = n_eqs)
    },
    
    #' @description Set options for EL evaluation.
    #' @template args-lambda_precision
    #' @template arg-support
    #' @template arg-lambda0
    set_opts = function(max_iter = 100, rel_tol = 1e-7, 
                        supp_corr = TRUE, lambda0 = 0) {
      private$.max_iter <- max_iter
      private$.rel_tol <- rel_tol
      private$.supp_corr <- supp_corr
      private$.lambda0 <- lambda0
    },
    
    #' @description Solve the inner optimization problem given G matrix.
    #' @template arg-G 
    #' @template arg-verbose
    lambda_nr = function(G, verbose = FALSE) {
      private$check_G(G)
      self$.G <- G
      lambdaNR(G = self$.G,
               lambda0 = private$.lambda0,
               support = private$.supp_corr,
               max_iter = private$.max_iter,
               rel_tol = private$.rel_tol,
               verbose = verbose)
    },
    
    #' @description Calculate log EL given G matrix.
    #' @template arg-G
    #' @param return_omega A boolean indicating whether to return the probability vector omega used to evaluate log EL.
    #' @param return_dldG A boolean indicating whether to return the gradient matrix dldG of log EL w.r.t. G.
    #' @template arg-verbose
    logel = function(G, return_omega = FALSE, return_dldG = FALSE, verbose = FALSE) {
      private$check_G(G)
      self$.G <- G
      logEL(G = self$.G, 
            lambda0 = private$.lambda0,
            max_iter = private$.max_iter,
            rel_tol = private$.rel_tol,
            support = private$.supp_corr,
            return_omega = return_omega, 
            return_dldG = return_dldG, 
            verbose = verbose)
    }
  )
)