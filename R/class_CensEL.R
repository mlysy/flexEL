#' @title Right-censored regression EL class
#' 
#' @description R6 class for EL regression with right-censored responses.
#' 
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
    .smooth = FALSE,
    .smooth_s = NULL,
    
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
    
    #' @field lambda0 Access or reset the initial value of lambda. The value 
    #'   must be a vector of length `n_eqs` (default to a vector of 0).
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
    
    #' @field max_iter_nr Access or reset the value of the maximum number of 
    #'   iterations for the Newton-Raphson algorithm. The value must be a 
    #'   positive integer (default to 100).
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
    
    #' @field rel_tol Access or reset the value of relative tolerance controlling 
    #'   the accuracy at convergence of the Newton-Raphson algorithm. The value 
    #'   must be a small positive number (default to 1e-7).
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
    
    #' @field max_iter_em Access or reset the value of the maximum number of 
    #'   iterations for the EM algorithm. The value must be a positive integer 
    #'   (default to 100).
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
    
    #' @field abs_tol Access or reset the value of absolute tolerance controlling 
    #'   the accuracy at convergence of the EM algorithm. The value must be 
    #'   a small positive number (default to 1e-3).
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
    
    #' @field supp_adj Access or reset the support correction flag. The value 
    #'   must be a boolean indicating whether to conduct support correction or 
    #'   not (default to FALSE).
    supp_adj = function(value) {
      if (missing(value)) private$.supp_adj
      else {
        value <- as.logical(value)
        if (is.na(value)) stop("Could not coerce `supp_adj` to boolean.")
        private$.supp_adj <- value
        private$.supp_adj_a <- max(1.0,0.5*log(CensEL_get_n_obs(private$.CEL)))
        CensEL_set_supp_adj(private$.CEL, private$.supp_adj, private$.supp_adj_a)
      }
    },
    
    #' @field supp_adj_a Access or reset the value of support corection factor.
    #'  The value must be a positive scalar (default to `max(1.0, log(n_obs)/2)`).
    supp_adj_a = function(value) {
      if (missing(value)) private$.supp_adj_a
      else if (!is.null(value) && (!is.numeric(value) || value <= 0)) {
        stop("`supp_adj_a` must be a positive number if not NULL.")
      }
      else {
        private$.supp_adj_a <- value
        CensEL_set_supp_adj(private$.CEL, private$.supp_adj, private$.supp_adj_a)
      }
    },
    
    #' @field smooth Access or reset the continuity correction flag. The value 
    #'   must be a boolean indicating whether to conduct continuity 
    #'   correction or not (default to FALSE). If set to `TRUE`, the tuning 
    #'   parameter for continuity correction is default to 10.
    smooth = function(value) {
      if (missing(value)) private$.smooth
      else {
        value <- as.logical(value)
        if (is.na(value)) stop("Could not coerce `smooth` to boolean.")
        private$.smooth <- value
        private$.smooth_s <- 10
        CensEL_set_smooth(private$.CEL, private$.smooth, private$.smooth_s)
      }
    },
    
    #' @field smooth_s Access or reset the continuity correction flag. The value 
    #'   must be a positive scalar for tuning the extent of 
    #'   continuity correction (default to 10). The smaller the value, the more 
    #'   smooth it makes.
    smooth_s = function(value) {
      if (missing(value)) private$.smooth_s
      else if (!is.null(value) && (!is.numeric(value) || value <= 0)) {
        stop("`smooth_s` must be a positive number if not NULL.")
      }
      else {
        private$.smooth_s <- value
        CensEL_set_smooth(private$.CEL, private$.smooth, private$.smooth_s)
      }
    }
  ),
  
  public = list(
    
    #' @description Create a new CensEL object.
    #' @param n_obs Number of observations.
    #' @param n_eqs Number of (moment constraint) equations.
    #' @return A `CensEL` object.
    initialize = function(n_obs, n_eqs) {
      if (n_obs <= 0 || n_eqs <= 0 || n_obs %% 1 != 0 || n_eqs %% 1 != 0) {
        stop("`n_obs` and `n_eqs` must be positive integers.")
      }
      private$.CEL <- CensEL_ctor(n_obs, n_eqs)
      private$.lambda0 <- rep(0, n_eqs)
    },
    
    #' @description Set the support correction flag and support correction factor.
    #' @param supp_adj     A boolean indicating whether to conduct support correction or not.
    #' @param supp_adj_a   Support adjustment factor. Defaults to `max(1.0, log(n_obs)/2)`.
    set_supp_adj = function(supp_adj = FALSE, supp_adj_a = NULL) {
      self$supp_adj <- supp_adj
      self$supp_adj_a <- supp_adj_a
    },
    
    #' @description Set the support correction flag and support correction factor.
    #' @param smooth     A boolean indicating whether to conduct support correction or not.
    #' @param smooth_s   Support adjustment factor. Defaults to `max(1.0, log(n_obs)/2)`.
    set_smooth = function(smooth = FALSE, smooth_s = NULL) {
      self$smooth <- smooth
      self$smooth_s <- smooth_s
    },
    
    #' @description Set more than one options together.
    #' @param max_iter_nr   A positive integer controlling the maximum number of iterations for the Newton-Raphson algorithm.
    #' @param rel_tol       A small positive number controlling accuracy at convergence for the Newton-Raphson algorithm.
    #' @param max_iter_em   A positive integer controlling the maximum number of iterations for the EM algorithm.
    #' @param abs_tol       A small positive number controlling accuracy at convergence for the EM algorithm.
    #' @param lambda0       Initial value of lambda of length `n_eqs`.
    #' @param supp_adj      A boolean indicating whether to conduct support correction or not.
    #' @param supp_adj_a    Support adjustment factor. Defaults to `max(1.0, log(n_obs)/2)`.
    #' @param smooth     A boolean indicating whether to conduct support correction or not.
    #' @param smooth_s   Support adjustment factor. Defaults to `max(1.0, log(n_obs)/2)`.
    set_opts = function(max_iter_nr = 100, rel_tol = 1e-7, 
                        max_iter_em = 100, abs_tol = 1e-3, 
                        supp_adj = FALSE, supp_adj_a = NULL,
                        smooth = FALSE, smooth_s = NULL, 
                        lambda0 = rep(0, CensEL_get_n_eqs(private$.CEL))) {
      self$max_iter_nr <- max_iter_nr
      self$rel_tol <- rel_tol
      self$max_iter_em <- max_iter_em
      self$abs_tol <- abs_tol
      self$set_supp_adj(supp_adj = supp_adj, supp_adj_a = supp_adj_a)
      self$set_smooth(smooth = smooth, smooth_s = smooth_s)
      self$lambda0 <- lambda0
    },
    
    #' @description Calculate the weights corresponding to the censored loglikelihood with EM algorithm.
    #' @param delta   A vector of censoring indicators of length `n_obs`, where 
    #'   0 means censored and 1 means observed. With support correction, the last 
    #'   entry of delta will be assumed to be 0.
    #' @param epsilon A vector of residuals of length `n_obs`. With support 
    #'   correction, the last entry of epsilon will be assumed to be -Inf.
    #' @param omega   A probability vector of length `n_obs + supp_adj`. with 
    #'   support correction, last entry of omega is assumed to correspond to 
    #'   the additional estimating equation.
    eval_weights = function(delta, epsilon, omega) {
      if (!all(delta %in% c(0,1))) {
        stop("`delta` must contain only 0 or 1s.")
      }
      if (sum(omega) != 1 && !all(omega >= 0 & omega <= 1)) {
        stop("`omega` must be a probability vector.")
      }
      if (!self$supp_adj && 
          (length(delta) != length(epsilon) || length(delta) != length(omega))) {
        stop("Without support correction, `delta`, `epsilon`, and `omega` must have the same length.")
      } else if (self$supp_adj &&
                 (length(delta) + 1 != length(omega) || length(epsilon) + 1 != length(omega))) {
        stop("With support correction, length of `delta` and `omega` must be the same, and 1 less than the length of `omega`.")
      }
      CensEL_eval_weights(private$.CEL, delta, epsilon, omega)
    },
    
    #' @description Calculate the probability vector base on the given G matrix, 
    #'   censoring indicator, and residuals.
    #' @param G       A matrix of dimension `n_obs x n_eqs`.
    #' @param delta   A vector of censoring indicators where 0 means censored and 1 means observed.
    #' @param epsilon A vector of residuals.
    omega_hat = function(G, delta, epsilon) {
      private$check_G(G)
      if (!all(delta %in% c(0,1))) {
        stop("`delta` must contain only 0 or 1s.")
      }
      if (length(delta) != nrow(G) || length(epsilon) != nrow(G)) {
        stop("`delta` and `epsilon` must have the same length and should be equal to the number of columns of `G`.")
      }
      CensEL_omega_hat(private$.CEL, t(G), delta, epsilon)
    },
    
    #' @description Calculate the log empirical likelihood base on the given G matrix.
    #' @param G       A matrix of dimension `n_obs x n_eqs`.
    #' @param delta   A vector of censoring indicators of length `n_obs`, where 
    #'   0 means censored and 1 means observed. With support correction, the last 
    #'   entry of delta will be assumed to be 0.
    #' @param epsilon A vector of residuals.
    #' @return A scalar.
    logel = function(G, delta, epsilon) {
      private$check_G(G)
      if (!all(delta %in% c(0,1))) {
        stop("`delta` must contain only 0 or 1s.")
      }
      if (length(delta) != nrow(G) || length(epsilon) != nrow(G) ) {
        stop("`delta` and `epsilon` must have the same length and should be equal to the number of columns of `G`.")
      }
      CensEL_logel(private$.CEL, t(G), delta, epsilon)
    }
  )
)