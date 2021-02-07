#' R6 class for EL with general moment specification.
#' 
#' A general EL object.
#' @export
GenEL <- R6::R6Class(
  
  classname = "GenEL",
  
  private = list(
    
    .GEL = NULL,
    .G = NULL,
    
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
    G = function(value) {
      if (missing(value)) private$.G
      else {
        private$check_G(value)
        private$.G <- value
      }
    }
  ),
  
  public = list(
    
    #' @description Create a new GenEL object.
    #' @param n_obs Number of observations.
    #' @param n_eqs Number of (moment constraint) equations.
    initialize = function(n_obs, n_eqs) {
      private$.GEL <- GenEL_ctor(n_obs, n_eqs)
      private$.G <- matrix(NA, nrow = n_obs, ncol = n_eqs)
    },
    
    set_max_iter = function(max_iter = 100) {
      if (!is.numeric(max_iter) | max_iter <= 0) {
        stop("`max_iter` must be a positive number.")
      }
      GenEL_set_max_iter(private$.GEL, max_iter)
    },
    
    set_rel_tol = function(rel_tol = 1e-7) {
      if (!is.numeric(rel_tol) | rel_tol <= 0) {
        stop("`rel_tol` must be a positive number.")
      }
      GenEL_set_rel_tol(private$.GEL, rel_tol)
    },
    
    set_supp_adj = function(supp_adj = FALSE, a = NULL) {
      if (!is.logical(supp_adj)) {
        stop("`supp_adj` must be a boolean.")
      }
      else if (!is.null(a)) {
        if (!is.numeric(a) | a <= 0) {
          stop("`a` must be a positive number if not NULL.")
        }
      }
      GenEL_set_supp_adj(private$.GEL, supp_adj, a)
    },
    
    set_lambda0 = function(lambda0) {
      # TODO: write a getter for n_obs and n_eqs in GenEL
      #   and check here that the length of lambda0 should be n_eqs (does supp_adj affects this?)
      GenEL_set_lambda0(private$.GEL, lambda0)
    }
    

    #' ,
    #' 
    #' #' @description Set options for EL evaluation.
    #' #' @template args-lambda_precision
    #' #' @template arg-support
    #' #' @template arg-lambda0
    #' set_opts = function(max_iter = 100, rel_tol = 1e-7, 
    #'                     supp_corr = TRUE, lambda0 = 0) {
    #'   private$.max_iter <- max_iter
    #'   private$.rel_tol <- rel_tol
    #'   private$.supp_corr <- supp_corr
    #'   private$.lambda0 <- lambda0
    #' },
    #' 
    #' #' @description Solve the inner optimization problem given G matrix.
    #' #' @template arg-G 
    #' #' @template arg-verbose
    #' lambda_nr = function(G, verbose = FALSE) {
    #'   private$check_G(G)
    #'   self$.G <- G
    #'   lambdaNR(G = self$.G,
    #'            lambda0 = private$.lambda0,
    #'            support = private$.supp_corr,
    #'            max_iter = private$.max_iter,
    #'            rel_tol = private$.rel_tol,
    #'            verbose = verbose)
    #' },
    #' 
    #' #' @description Calculate log EL given G matrix.
    #' #' @template arg-G
    #' #' @param return_omega A boolean indicating whether to return the probability vector omega used to evaluate log EL.
    #' #' @param return_dldG A boolean indicating whether to return the gradient matrix dldG of log EL w.r.t. G.
    #' #' @template arg-verbose
    #' logel = function(G, return_omega = FALSE, return_dldG = FALSE, verbose = FALSE) {
    #'   private$check_G(G)
    #'   self$.G <- G
    #'   logEL(G = self$.G, 
    #'         lambda0 = private$.lambda0,
    #'         max_iter = private$.max_iter,
    #'         rel_tol = private$.rel_tol,
    #'         support = private$.supp_corr,
    #'         return_omega = return_omega, 
    #'         return_dldG = return_dldG, 
    #'         verbose = verbose)
    #' }
  )
)