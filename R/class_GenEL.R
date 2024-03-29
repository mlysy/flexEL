#' @title General EL class
#'
#' @description R6 class for EL with general moment specification.
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
    .weight_adj = NULL,

    #' @description Check the dimension of G is what is expected.
    #' @template arg-G
    check_G = function(G) {
      n_obs <- GenEL_get_n_obs(private$.GEL)
      n_eqs <- GenEL_get_n_eqs(private$.GEL)
      if (nrow(G) != n_obs | ncol(G) != n_eqs) {
        stop("Dimension of `G` matrix must be `n_obs` x `n_eqs`.")
      }
    },

    #' @description Check the dimension of weights is as expected.
    #' @param weights A numeric vector.
    check_weights = function(weights) {
      if (length(weights) != GenEL_get_n_obs(private$.GEL)) {
        stop("Length of `weights` does not equal to the number of obserations.")
      }
      if (any(weights < 0)) {
        stop("`weights` should contain only non-negative values.")
      }
    }
  ),

  active = list(

    #' @field max_iter Access or reset the value of the maximum number of
    #'   iterations. The value must be a positive integer (default to 100).
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

    #' @field rel_tol Access or reset the value of relative tolerance
    #'   controlling the accuracy at convergence. The value must be a small
    #'   positive number (default to 1e-7).
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

    #' @field lambda0 Access or reset the initial value of lambda. The value
    #'   must be a vector of length `n_eqs` (default to a vector of 0).
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

    #' @field supp_adj Access or reset the support correction flag. The value
    #'   must be a boolean indicating whether to conduct support correction or
    #'   not (default to FALSE).
    supp_adj = function(value) {
      if (missing(value)) private$.supp_adj
      else {
        value <- as.logical(value)
        if(is.na(value)) stop("Could not coerce `supp_adj` to boolean.")
        private$.supp_adj <- value
        GenEL_set_supp_adj(private$.GEL,
                           private$.supp_adj,
                           private$.supp_adj_a,
                           private$.weight_adj)
      }
    },

    #' @field supp_adj_a Access or reset the value of support corection factor.
    #'   The value must be a scalar (defaults to `max(1.0, log(n_obs)/2)`).
    supp_adj_a = function(value) {
      if (missing(value)) private$.supp_adj_a
      else if (!is.null(value) && (!is.numeric(value) || value <= 0)) {
        stop("`supp_adj_a` must be a positive number if not NULL.")
      }
      else {
        private$.supp_adj_a <- value
        GenEL_set_supp_adj(private$.GEL,
                           private$.supp_adj,
                           private$.supp_adj_a,
                           private$.weight_adj)
      }
    },

    #' @field weight_adj Access or reset the value of the weight for the
    #'   additional fake observation under support correction. The value must
    #'   be a scalar (default to NULL).
    weight_adj = function(value) {
      if (missing(value)) private$.weight_adj
      else if (!is.null(value) && (!is.numeric(value) || value <= 0)) {
        stop("`weight_adj` must be a positive number if not NULL.")
      }
      else {
        private$.weight_adj <- value
        GenEL_set_supp_adj(private$.GEL,
                           private$.supp_adj,
                           private$.supp_adj_a,
                           private$.weight_adj)
      }
    }
  ),

  public = list(

    #' @description Create a new `GenEL` object.
    #' @param n_obs Number of observations.
    #' @param n_eqs Number of (moment constraint) equations.
    #' @return A `GenEL` object.
    initialize = function(n_obs, n_eqs) {
      private$.GEL <- GenEL_ctor(n_obs, n_eqs)
      private$.lambda0 <- rep(0, n_eqs)
      private$.supp_adj_a <- max(1.0,0.5*log(n_obs))
      private$.weight_adj <- 1.0
    },

    #' @description Set the support correction flag and support correction factor.
    #' @param supp_adj     A boolean indicating whether to conduct support correction or not.
    #' @param supp_adj_a   Support adjustment factor (default to `max(1.0, log(n_obs)/2)`).
    #' @param weight_adj   Weight under weighted log EL (default to 1.0).
    set_supp_adj = function(supp_adj = FALSE, supp_adj_a = NULL, weight_adj = NULL) {
      self$supp_adj <- supp_adj
      self$supp_adj_a <- supp_adj_a
      self$weight_adj <- weight_adj
    },

    #' @description Set more than one options together.
    #' @param max_iter   A positive integer controlling the maximum number of iterations.
    #' @param rel_tol    A small positive number controlling accuracy at convergence.
    #' @param supp_adj   A boolean indicating whether to conduct support correction or not.
    #' @param supp_adj_a Support adjustment factor. Defaults to `max(1.0, log(n_obs)/2)`.
    #' @param weight_adj   Weight under weighted log EL (default to 1.0).
    #' @param lambda0    Initialization vector of size `n_eqs`.  Defaults to a vector of zeros.
    set_opts = function(max_iter = 100, rel_tol = 1e-7,
                        supp_adj = FALSE, supp_adj_a = NULL, weight_adj = NULL,
                        lambda0 = NULL) {
      self$max_iter <- max_iter
      self$rel_tol <- rel_tol
      self$set_supp_adj(supp_adj = supp_adj, supp_adj_a = supp_adj_a, weight_adj = weight_adj)
      if(is.null(lambda0)) lambda0 <- rep(0, GenEL_get_n_eqs(private$.GEL))
      self$lambda0 <- lambda0
    },

    #' @description Calculate the solution of the dual problem of maximum log EL problem.
    #' @param G        A numeric matrix of dimension `n_obs x n_eqs`.
    #' @param weights  A numeric vector of length `n_obs` containing non-negative values.
    #' @param check_conv  If `TRUE`, checks whether the desired `rel_tol` was reached within the given maximum number of Newton-Raphson iterations.  If not, returns a vector of `NaN`s.
    #' @return A numeric vector of length `n_eqs`.
    lambda_nr = function(G, weights, check_conv = TRUE) {
      private$check_G(G)
      n_obs <- GenEL_get_n_obs(private$.GEL)
      if (missing(weights) || is.null(weights)) {
        weights <- rep(1.0, n_obs)
      }
      private$check_weights(weights)
      GenEL_lambda_nr(private$.GEL, t(G), weights, check_conv)
    },

    #' @description Calculate the probability vector base on the given G matrix.
    #' @param G        A numeric matrix of dimension `n_obs x n_eqs`.
    #' @param weights  A numeric vector of length `n_obs` containing non-negative values.
    #' @param check_conv  If `TRUE`, checks whether the desired `rel_tol` was reached within the given maximum number of Newton-Raphson iterations.  If not, returns a vector of `NaN`s.
    #' @return A probability vector of length `n_obs + supp_adj`.
    omega_hat = function(G, weights, check_conv = TRUE) {
      private$check_G(G)
      n_obs <- GenEL_get_n_obs(private$.GEL)
      if (missing(weights) || is.null(weights)) {
        weights <- rep(1.0, n_obs)
      }
      ## if (length(weights) != n_obs) {
      ##   stop("Length of `weights` does not equal to the number of obserations.")
      ## }
      ## if (any(weights < 0)) {
      ##   stop("`weights` should contain only non-negative values.")
      ## }
      private$check_weights(weights)
      lambda <- self$lambda_nr(G = G, weights = weights,
                               check_conv = check_conv)
      GenEL_omega_hat(private$.GEL, lambda, t(G), weights)
    },

    #' @description Calculate the log empirical likelihood for a given `G` matrix.
    #' @param G       A numeric matrix of dimension `n_obs x n_eqs`.
    #' @param weights  A numeric vector of length `n_obs` containing non-negative values.
    #' @param check_conv  If `TRUE`, checks whether the desired `rel_tol` was reached within the given maximum number of Newton-Raphson iterations.  If not, returns `-Inf`.
    #' @return A scalar.
    logel = function(G, weights, check_conv = TRUE) {
      private$check_G(G)
      if (missing(weights) || is.null(weights)) {
        ans <- GenEL_logel(private$.GEL, t(G), check_conv)
      }
      else {
        ## n_obs <- GenEL_get_n_obs(private$.GEL)
        ## if (length(weights) != n_obs) {
        ##   stop("Length of `weights` does not equal to the number of obserations.")
        ## }
        ## if (any(weights < 0)) {
        ##   stop("`weights` should contain only non-negative values.")
        ## }
        private$check_weights(weights)
        ans <- GenEL_weighted_logel(private$.GEL, t(G), weights, check_conv)
      }
      ans
    },

    #' @description Calculate log EL and the derivative of log EL w.r.t. G evaluated at G.
    #' @param G        A numeric matrix of dimension `n_obs x n_eqs`.
    #' @param weights  A numeric vector of length `n_obs` containing non-negative values.
    #' @param check_conv If `TRUE`, checks whether the desired `rel_tol` was reached within the given maximum number of Newton-Raphson iterations.  If not, sets the log EL to negative infinity and the gradient to a matrix of `NaN`s.
    #' @return A list with elements `logel` and `grad`, the latter being a matrix of size `n_obs x n_eqs`.
    logel_grad = function(G, weights, check_conv = TRUE) {
      private$check_G(G)
      if (missing(weights) || is.null(weights)) {
        ans <- GenEL_logel_grad(private$.GEL, t(G), check_conv)
      }
      else {
        ## n_obs <- GenEL_get_n_obs(private$.GEL)
        ## if (length(weights) != n_obs) {
        ##   stop("Length of `weights` does not equal to the number of obserations.")
        ## }
        ## if (any(weights < 0)) {
        ##   stop("`weights` should contain only non-negative values.")
        ## }
        private$check_weights(weights)
        ans <- GenEL_weighted_logel_grad(private$.GEL, t(G), weights, check_conv)
      }
      ans$grad <- t(ans$grad)
      ans
    }
  )
)
