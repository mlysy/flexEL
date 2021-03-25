#' Constructor and methods for general EL objects.
#' 
#' R6 class for EL with general moment specification.
#' 
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
    
    #' @description Access or reset the value of the maximum number of iterations.
    #' @param value Missing or a positive integer (default to 100).
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
    
    #' @description Access or reset the value of relative tolerance controlling the accuracy at convergence..
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
    
    #' @description Access or reset the initial value of lambda.
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
    
    #' @description Access or reset the support correction flag.
    #' @param value Missing or a boolean indicating whether to conduct support correction or not.
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
    
    #' @description Access or reset the value of support corection factor.
    #' @param value Missing or a scalar. Defaults to `max(1.0, log(n_obs)/2)`.
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
    #' @return A `GenEL` object.
    initialize = function(n_obs, n_eqs) {
      private$.GEL <- GenEL_ctor(n_obs, n_eqs)
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
      GenEL_set_supp_adj(private$.GEL, supp_adj, supp_adj_a)
    },
    
    #' @description Set more than one options together.
    #' @param max_iter   A positive integer controlling the maximum number of iterations.
    #' @param rel_tol    A small positive number controlling accuracy at convergence.
    #' @param lambda0    Initialization vector of size `n_eqs`.
    #' @param supp_adj   A boolean indicating whether to conduct support correction or not.
    #' @param supp_adj_a Support adjustment factor. Defaults to `max(1.0, log(n_obs)/2)`.
    set_opts = function(max_iter = 100, rel_tol = 1e-7, 
                        lambda0 = rep(0, GenEL_get_n_eqs(private$.GEL)), 
                        supp_adj = FALSE, supp_adj_a = NULL) {
      self$max_iter <- max_iter
      self$rel_tol <- rel_tol
      self$lambda0 <- lambda0
      self$set_supp_adj(supp_adj = supp_adj, supp_adj_a = supp_adj_a)
    },
    
    #' @description Calculate the solution of the dual problem of maximum log EL problem.
    #' @param G        A matrix of dimension `n_eqs x n_obs`.
    #' @param verbose  A boolean indicating whether to print out number of iterations and maximum error at the end of the Newton-Raphson algorithm.
    #' @return A numeric vector of length `n_eqs`.
    lambda_nr = function(G, verbose = FALSE) {
      private$check_G(G)
      GenEL_lambda_nr(private$.GEL, t(G), verbose)
    },
    
    #' @description Calculate the probability vector base on the given G matrix.
    #' @param G        A matrix of dimension `n_obs x n_eqs`.
    #' @param verbose  A boolean indicating whether to print out number of iterations and maximum error at the end of the Newton-Raphson algorithm.
    #' @return A probability vector of length `n_obs + supp_adj`.
    omega_hat = function(G, verbose = FALSE) {
      private$check_G(G)
      lambda <- self$lambda_nr(G, verbose)
      GenEL_omega_hat(private$.GEL, lambda, t(G))
    },
    
    #' @description Calculate the log empirical likelihood base on the given probability vector omega.
    #' @param omega       A probability vector of length `n_eqs + supp_adj`.
    #' @return A scalar.
    logel_omega = function(omega) {
      n_obs2 <- GenEL_get_n_obs(private$.GEL) + GenEL_supp_adj(private$.GEL)
      if (length(omega) != n_obs2) {
        stop("`omega` must have length `n_eqs + supp_adj`.")
      }
      GenEL_logel_omega(private$.GEL, omega)
    },
    
    #' @description Calculate the log empirical likelihood base on the given G matrix.
    #' @param G       A matrix of dimension `n_eqs x n_obs`.
    #' @return A scalar.
    logel = function(G) {
      GenEL_logel(G)
    },
    
    #' @description Calculate the probability vector, log EL, and the derivative of log EL w.r.t. G evaluated at G.
    #' @param G        A matrix of dimension `n_eqs x n_obs`.
    #' @param verbose  A boolean indicating whether to print out number of iterations and maximum error at the end of the Newton-Raphson algorithm.
    #' @return A list of three elements.
    logel_grad = function(G, verbose = FALSE) {
      private$check_G(G)
      GenEL_Logel_grad(private$.GEL, G, verbose)
    }
  )
)