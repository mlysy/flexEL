#' R6 class for EL with general moment specification.
#' 
#' A general EL object.
#' @export
GenEL <- R6::R6Class(
  
  classname = "GenEL",
  
  private = list(
    
    .GEL = NULL,
    .lambda0 = NULL,
    .max_iter = 100,
    .rel_tol = 1e-7,
    .supp_adj = FALSE,
    .supp_adj_a = NULL,
    
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
    
    #' @description Access or reset the value of .max_iter.
    #' @param value Missing or a positive integer (default to 200).
    max_iter = function(value) {
      if (missing(value)) private$.max_iter
      else if (!is.numeric(value) | value <= 0) {
        stop("`max_iter` must be a positive number.")
      }
      else {
        private$.max_iter <- value
        GenEL_set_max_iter(private$.GEL, private$.max_iter)
      }
    },
    
    #' @description Access or reset the value of .rel_tol.
    #' @param value Missing or a small positive number (default to 1e-7).
    rel_tol = function(value) {
      if (missing(value)) private$.rel_tol
      else if (!is.numeric(value) | value <= 0) {
        stop("`rel_tol` must be a positive number.")
      }
      else {
        private$.rel_tol <- value
        GenEL_set_rel_tol(private$.GEL, private$.rel_tol)
      }
    },
    
    #' @description Access or reset the value of .lambda0.
    #' @param value Missing or a vector of length `n_eqs`.
    lambda0 = function(value) {
      if (missing(value)) private$.lambda0
      else {
        n_eqs <- GenEL_get_n_eqs(private$.GEL)
        if (length(value) != n_eqs) {
          stop("`lambda0` must be of length equal to the number of estimating equations.")
        }
        private$.lambda0 <- value
        GenEL_set_lambda0(private$.GEL, private$.lambda0)
      }
    },
    
    supp_adj = function(value) {
      if (missing(value)) private$.supp_adj
      else if (!is.logical(supp_adj)) {
        stop("`supp_adj` must be a boolean.")
      }
      else {
        private$.supp_adj <- value
        private$.supp_adj_a <- max(1.0,0.5*log(GenEL_get_n_obs(private$.GEL)))
        GenEL_set_supp_adj(private$.GEL, private$.supp_adj, private$.supp_adj_a)
      }
    },
    
    supp_adj_a = function(value) {
      if (missing(value)) private$.supp_adj_a
      else if (!is.numeric(value) | value <= 0) {
        stop("`a` must be a positive number if not NULL.")
      }
      else {
        private$.supp_adj_a <- value
        GenEL_set_supp_adj(private$.GEL, private$.supp_adj, private$.supp_adj_a)
      }
    }
  ),
  
  public = list(
    
    #' @description Create a new GenEL object.
    #' @param n_obs Number of observations.
    #' @param n_eqs Number of (moment constraint) equations.
    initialize = function(n_obs, n_eqs) {
      private$.GEL <- GenEL_ctor(n_obs, n_eqs)
      private$.lambda0 <- rep(0, n_eqs)
    },
    
    set_supp_adj = function(supp_adj = FALSE, supp_adj_a = NULL) {
      if (!is.logical(supp_adj)) {
        stop("`supp_adj` must be a boolean.")
      }
      else if (!is.null(supp_adj_a)) {
        if (!is.numeric(supp_adj_a) | supp_adj_a <= 0) {
          stop("`a` must be a positive number if not NULL.")
        }
      }
      GenEL_set_supp_adj(private$.GEL, supp_adj, supp_adj_a)
    },
    
    set_opts = function(max_iter = 100, rel_tol = 1e-7, 
                        lambda0 = rep(0, GenEL_get_n_eqs(private$.GEL)), 
                        supp_adj = FALSE, a = NULL) {
      self$max_iter <- max_iter
      self$rel_tol <- rel_tol
      self$lambda0 <- lambda0
      self$set_supp_adj(supp_adj = supp_adj, a = a)
    },
    
    lambda_nr = function(G, verbose = FALSE) {
      private$check_G(G)
      GenEL_lambda_nr(private$.GEL, G, verbose)
    },
    
    omega_hat = function(G, verbose = FALSE) {
      private$check_G(G)
      lambda <- self$lambda_nr(G, verbose)
      GenEL_omega_hat(private$.GEL, lambda, G)
    },
    
    logel = function(G, verbose = FALSE) {
      private$check_G(G)
      omega <- self$omega_hat(G, verbose)
      GenEL_logel_omega(private$.GEL, omega)
    },
    
    logel_grad = function(G, verbose = FALSE) {
      private$check_G(G)
      GenEL_Logel_grad(private$.GEL, G, verbose)
    }
  )
)