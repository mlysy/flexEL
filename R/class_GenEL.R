#' R6 class for EL with general moment specification.
#' 
#' A general EL object.
#' @export
GenEL <- R6::R6Class(
  
  classname = "GenEL",
  
  private = list(
    
    .GEL = NULL,
    .G = NULL,
    .lambda = NULL,
    .omega = NULL,
    
    #' @description Check the dimension of G is what is expected.
    #' @param G A numeric matrix.
    check_G = function(G) {
      n_obs <- GenEL_get_n_obs(private$.GEL)
      n_eqs <- GenEL_get_n_eqs(private$.GEL)
      if (nrow(G) != n_obs | ncol(G) != n_eqs) {
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
      n_eqs <- GenEL_get_n_eqs(private$.GEL)
      if (length(lambda0) != n_eqs) {
        stop("`lambda0` must be of length equal to the number of estimating equations.")
      }
      GenEL_set_lambda0(private$.GEL, lambda0)
    },
    
    lambda_nr = function() {
      if (is.null(private$.G)) {
        stop("G matrix has not been set.")
      }
      .lambda <- GenEL_lambda_nr(private$.GEL, private$.G)
    },
    
    omega_hat = function() {
      if (is.null(private$.G)) {
        stop("G matrix has not been set.")
      }
      if (is.null(private$.lambda)) {
        stop("lamba has not been calculated.")
      }
      .omega <- GenEL_omega_hat(private$.GEL, private$.lambda, private$.G)
    },
    
    logel = function() {
      if (is.null(private$.omega)) {
        stop("omega has not been calculated.")
      }
      GenEL_logel_omega(private$.GEL, private$.omega)
    }
  )
)