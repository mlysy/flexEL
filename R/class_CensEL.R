#' R6 class for EL regression where response is under right censoring.
#' 
#' A general EL object.
#' @export
CensEL <- R6::R6Class(
  
  classname = "CensEL",
  
  private = list(
    .CEL = NULL,
    .lambda0 = NULL,
    .max_iter_nr = 100,
    .rel_tol = 1e-7,
    .max_iter_em = 100,
    .abs_tol = 1e-3,
    .supp_adj = FALSE,
    .supp_adj_a = NULL,
    
    #' @description Check the dimension of G is what is expected.
    #' @param G A numeric matrix.
    check_G = function(G) {
      n_obs <- CensEL_get_n_obs(private$.CEL)
      n_eqs <- CensEL_get_n_eqs(private$.CEL)
      if (nrow(G) != n_obs | ncol(G) != n_eqs) {
        stop("Dimension of G matrix must be `n_obs` x `n_eqs`.")
      }
    }
  ),
  
  active = list(
    
    #' @description Access or reset the initial value of lambda.
    #' @param value Missing or a vector of length `n_eqs`.
    lambda0 = function(value) {
      if (missing(value)) private$.lambda0
      else {
        n_eqs <- CensEL_get_n_eqs(private$.CEL)
        if (length(value) != n_eqs) {
          stop("`lambda0` must be of length equal to the number of estimating equations.")
        }
        private$.lambda0 <- value
        CensEL_set_lambda0(private$.CEL, private$.lambda0)
      }
    },
    
    #' @description Access or reset the value of the maximum number of iterations for the Newton-Raphson algorithm.
    #' @param value Missing or a positive integer (default to 100).
    max_iter_nr = function(value) {
      if (missing(value)) private$.max_iter_nr
      else if (!is.numeric(value) | value <= 0) {
        stop("`max_iter` must be a positive number.")
      }
      else {
        private$.max_iter_nr <- value
        CensEL_set_max_iter_nr(private$.CEL, private$.max_iter_nr)
      }
    },
    
    #' @description Access or reset the value of relative tolerance controlling the accuracy at convergence of the Newton-Raphson algorithm.
    #' @param value Missing or a small positive number (default to 1e-7).
    rel_tol = function(value) {
      if (missing(value)) private$.rel_tol
      else if (!is.numeric(value) | value <= 0) {
        stop("`rel_tol` must be a positive number.")
      }
      else {
        private$.rel_tol <- value
        CensEL_set_rel_tol(private$.CEL, private$.rel_tol)
      }
    },
    
    #' @description Access or reset the value of the maximum number of iterations for the EM algorithm.
    #' @param value Missing or a positive integer (default to 100).
    max_iter_em = function(value) {
      if (missing(value)) private$.max_iter_em
      else if (!is.numeric(value) | value <= 0) {
        stop("`max_iter` must be a positive number.")
      }
      else {
        private$.max_iter_em <- value
        CensEL_set_max_iter_em(private$.CEL, private$.max_iter_em)
      }
    },
    
    #' @description Access or reset the value of absolute tolerance controlling the accuracy at convergence of the EM algorithm.
    #' @param value Missing or a small positive number (default to 1e-3).
    abs_tol = function(value) {
      if (missing(value)) private$.abs_tol
      else if (!is.numeric(value) | value <= 0) {
        stop("`abs_tol` must be a positive number.")
      }
      else {
        private$.abs_tol <- value
        CensEL_set_abs_tol(private$.CEL, private$.abs_tol)
      }
    },
    
    #' @description Access or reset the support correction flag.
    #' @param value Missing or a boolean indicating whether to conduct support correction or not.
    supp_adj = function(value) {
      if (missing(value)) private$.supp_adj
      else if (!is.logical(value)) {
        stop("`supp_adj` must be a boolean.")
      }
      else {
        private$.supp_adj <- value
        private$.supp_adj_a <- max(1.0,0.5*log(CensEL_get_n_obs(private$.CEL)))
        CensEL_set_supp_adj(private$.CEL, private$.supp_adj, private$.supp_adj_a)
      }
    },
    
    #' @description Access or reset the value of support corection factor.
    #' @param value Missing or a scalar. Defaults to `max(1.0, log(n_obs)/2)`.
    supp_adj_a = function(value) {
      if (missing(value)) private$.supp_adj_a
      else if (!is.numeric(value) | value <= 0) {
        stop("`a` must be a positive number if not NULL.")
      }
      else {
        private$.supp_adj_a <- value
        CensEL_set_supp_adj(private$.CEL, private$.supp_adj, private$.supp_adj_a)
      }
    }
  ),
  
  public = list(
    
    #' @description Create a new CensEL object.
    #' @param n_obs Number of observations.
    #' @param n_eqs Number of (moment constraint) equations.
    #' @return A `CensEL` object.
    initialize = function(n_obs, n_eqs) {
      private$.CEL <- CensEL_ctor(n_obs, n_eqs)
      private$.lambda0 <- rep(0, n_eqs)
    },
    
    #' @description Set the support correction flag and support correction factor.
    #' @param supp_adj     A boolean indicating whether to conduct support correction or not.
    #' @param supp_adj_a   Support adjustment factor. Defaults to `max(1.0, log(n_obs)/2)`.
    set_supp_adj = function(supp_adj = FALSE, supp_adj_a = NULL) {
      if (!is.logical(supp_adj)) {
        stop("`supp_adj` must be a boolean.")
      }
      else if (!is.null(supp_adj_a)) {
        if (!is.numeric(supp_adj_a) | supp_adj_a <= 0) {
          stop("`a` must be a positive number if not NULL.")
        }
      }
      CensEL_set_supp_adj(private$.CEL, supp_adj, supp_adj_a)
    },
    
    #' @description Set more than one options together.
    #' @param max_iter_nr   A positive integer controlling the maximum number of iterations for the Newton-Raphson algorithm.
    #' @param rel_tol       A small positive number controlling accuracy at convergence for the Newton-Raphson algorithm.
    #' @param max_iter_em   A positive integer controlling the maximum number of iterations for the EM algorithm.
    #' @param abs_tol       A small positive number controlling accuracy at convergence for the EM algorithm.
    #' @param lambda0       Initialization vector of size `n_eqs`.
    #' @param supp_adj      A boolean indicating whether to conduct support correction or not.
    #' @param supp_adj_a    Support adjustment factor. Defaults to `max(1.0, log(n_obs)/2)`.
    set_opts = function(max_iter_nr = 100, rel_tol = 1e-7, 
                        max_iter_em = 100, abs_tol = 1e-3, 
                        lambda0 = rep(0, CensEL_get_n_eqs(private$.CEL)), 
                        supp_adj = FALSE, supp_adj_a = NULL) {
      self$max_iter_nr <- max_iter_nr
      self$rel_tol <- rel_tol
      self$max_iter_em <- max_iter_em
      self$abs_tol <- abs_tol
      self$lambda0 <- lambda0
      self$set_supp_adj(supp_adj = supp_adj, supp_adj_a = supp_adj_a)
    },
    
    eval_weights = function(delta, epsilon, omega) {
      CensEL_eval_weights(private$.CEL, delta, epsilon, omega)
    },
    
    omega_hat = function(G, delta, epsilon) {
      private$check_G(G)
      CensEL_omega_hat(private$.CEL, t(G), delta, epsilon)
    },
    
    #' @description Calculate the log empirical likelihood base on the given G matrix.
    #' @param G       A matrix of dimension `n_eqs x n_obs`.
    #' @return A scalar.
    logel = function(G, delta, epsilon) {
      private$check_G(G)
      CensEL_logel(private$.CEL, t(G), delta, epsilon)
    }
    
  )
  
)