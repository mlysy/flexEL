#' R6 Class for general EL models.
#'
#' @export
GenEL <- R6::R6Class(

  classname = "GenEL",

  private = list(

    .gelptr = NULL,
    .lambda0 = NULL,
    .max_iter = NULL,
    .rel_tol = NULL,
    .supp_adj = NULL,
    .supp_adj_a = NULL,
    .weight_adj = NULL,

    #' @description Check the dimension of G is what is expected.
    #' @template param_G
    check_G = function(G) {
      n_obs <- GenEL_get_n_obs(private$.gelptr)
      n_eqs <- GenEL_get_n_eqs(private$.gelptr)
      if(!is.numeric(G) || (nrow(G) != n_obs) || (ncol(G) != n_eqs)) {
        stop("`G` must be a numeric matrix of size `n_obs x n_eqs`.")
      }
    },

    #' @description Check the dimension of weights is as expected.
    #' @param weights A numeric vector.
    check_weights = function(weights) {
      if (length(weights) != GenEL_get_n_obs(private$.gelptr)) {
        stop("Length of `weights` must equal `n_obs`.")
      }
      if (!is.numeric(weights) || any(weights < 0)) {
        stop("`weights` must contain only non-negative values.")
      }
    },

    #' @description Get `private$.gelptr`.
    #' @details This is used to treat `.gelptr` as a protected member, i.e., for classes inheriting from `GenEL` to pass `.gelptr` on to C++ code.  Accessor is private as the intent is for this to not be called directly by the user.
    #' @note The returned `externalptr` is pointing to the same C++ object as is `private$.gelptr$.  So, modifying one modifies the other.
    get_gelptr = function() {
      private$.gelptr
    }
  ),

  active = list(

    #' @field n_obs Number of observations.
    n_obs = function() {
      GenEL_get_n_obs(private$.gelptr)
    },

    #' @field n_eqs Number of observations.
    n_eqs = function() {
      GenEL_get_n_eqs(private$.gelptr)
    },

    #' @field max_iter Maximum number of Newton-Raphson iterations.
    max_iter = function(value) {
      if (missing(value)) private$.max_iter
      else if (!is.numeric(value) || value <= 0) {
        stop("`max_iter` must be a positive number.")
      }
      else {
        private$.max_iter <- value
        GenEL_set_max_iter(private$.gelptr, as.integer(private$.max_iter))
      }
    },

    #' @field rel_tol Relative tolerance for Newton-Raphson algorithm.
    rel_tol = function(value) {
      if (missing(value)) {
        private$.rel_tol
      } else if (!is.numeric(value) || value <= 0) {
        stop("`rel_tol` must be a positive number.")
      } else {
        private$.rel_tol <- value
        GenEL_set_rel_tol(private$.gelptr, private$.rel_tol)
      }
    },

    #' @field lambda0 Initialization vector for the Newton-Raphson algorithm.
    lambda0 = function(value) {
      if (missing(value)) private$.lambda0
      else {
        n_eqs <- GenEL_get_n_eqs(private$.gelptr)
        if ((!is.numeric(value)) || (length(value) != n_eqs)) {
          stop("`lambda0` must be a numeric vector of length `n_eqs`.")
        }
        private$.lambda0 <- value
        GenEL_set_lambda0(private$.gelptr, private$.lambda0)
      }
    },

    #' @field supp_adj Whether or not to perform support adjustment.
    supp_adj = function(value) {
      if (missing(value)) private$.supp_adj
      else {
        value <- as.logical(value)
        if(is.na(value)) stop("Could not coerce `supp_adj` to boolean.")
        private$.supp_adj <- value
        GenEL_set_supp_adj(private$.gelptr,
                           private$.supp_adj,
                           private$.supp_adj_a,
                           private$.weight_adj)
      }
    },

    #' @field supp_adj_a Support adjustment factor.
    supp_adj_a = function(value) {
      if (missing(value)) private$.supp_adj_a
      else if (!is.null(value) && (!is.numeric(value) || value <= 0)) {
        stop("`supp_adj_a` must be a positive number if not NULL.")
      }
      else {
        private$.supp_adj_a <- value
        GenEL_set_supp_adj(private$.gelptr,
                           private$.supp_adj,
                           private$.supp_adj_a,
                           private$.weight_adj)
      }
    },

    #' @field weight_adj Weights of pseudo-observation for support adjustment.
    weight_adj = function(value) {
      if (missing(value)) private$.weight_adj
      else if (!is.numeric(value) || value <= 0) {
        stop("`weight_adj` must be a positive scalar.")
      }
      else {
        private$.weight_adj <- value
        GenEL_set_supp_adj(private$.gelptr,
                           private$.supp_adj,
                           private$.supp_adj_a,
                           private$.weight_adj)
      }
    }
  ),

  public = list(

    #' @description Create a new `GenEL` object.
    #'
    #' @template param_n_obs
    #' @template param_n_eqs
    #' @template param_max_iter_nr
    #' @template param_rel_tol
    #' @template param_supp_adj
    #' @template param_supp_adj_a
    #' @template param_weight_adj
    #' @template param_lambda0
    #'
    #' @return A `GenEL` object.
    initialize = function(n_obs, n_eqs,
                          max_iter = 100,
                          rel_tol = 1e-7,
                          supp_adj = FALSE,
                          supp_adj_a = max(1.0, .5 * log(n_obs)),
                          weight_adj = 1.0,
                          lambda0 = rep(0, n_eqs)) {
      private$.gelptr <- GenEL_ctor(n_obs, n_eqs)
      self$set_opts(max_iter = max_iter, rel_tol = rel_tol,
                    supp_adj = supp_adj, supp_adj_a = supp_adj_a,
                    weight_adj = weight_adj, lambda0 = lambda0)
    },

    ## #' @description Set the support correction flag and support correction factor.
    ## #' @param supp_adj     A boolean indicating whether to conduct support correction or not.
    ## #' @param supp_adj_a   Support adjustment factor (default to `max(1.0, log(n_obs)/2)`).
    ## #' @param weight_adj   Weight under weighted log EL (default to 1.0).
    ## set_supp_adj = function(supp_adj = FALSE, supp_adj_a = NULL, weight_adj = NULL) {
    ##   self$supp_adj <- supp_adj
    ##   self$supp_adj_a <- supp_adj_a
    ##   self$weight_adj <- weight_adj
    ## },

    #' @description Set multiple options together.
    #'
    #' @template param_max_iter_nr
    #' @template param_rel_tol
    #' @template param_supp_adj
    #' @template param_supp_adj_a
    #' @template param_weight_adj
    #' @template param_lambda0
    #'
    #' @details Any number of these arguments can be supplied.  Those that are missing are not modified.
    set_opts = function(max_iter, rel_tol,
                        supp_adj, supp_adj_a, weight_adj,
                        lambda0) {
      if(!missing(max_iter)) {
        self$max_iter <- max_iter
      }
      if(!missing(rel_tol)) {
        self$rel_tol <- rel_tol
      }
      if(!missing(supp_adj)) {
        self$supp_adj <- supp_adj
      }
      if(!missing(supp_adj_a)) {
        self$supp_adj_a <- supp_adj_a
      }
      if(!missing(weight_adj)) {
        self$weight_adj <- weight_adj
      }
      if(!missing(lambda0)) {
        self$lambda0 <- lambda0
      }
    },

    #' @description Calculate the solution of the empirical likelihood dual problem using the Newton-Raphson algorithm.
    #'
    #' @template param_G
    #' @template param_weights
    #' @template param_check_conv
    #'
    #' @return A numeric vector of length `n_eqs`.
    lambda_nr = function(G, weights, check_conv = TRUE) {
      private$check_G(G)
      n_obs <- GenEL_get_n_obs(private$.gelptr)
      if (missing(weights) || is.null(weights)) {
        weights <- rep(1.0, n_obs)
      }
      private$check_weights(weights)
      GenEL_lambda_nr(private$.gelptr, t(G), weights, check_conv)
    },

    #' @description Calculate the probability vector base on the given G matrix.
    #'
    #' @template param_G
    #' @template param_weights
    #' @template param_check_conv
    #'
    #' @return A probability vector of length `n_obs + supp_adj`.
    omega_hat = function(G, weights, check_conv = TRUE) {
      private$check_G(G)
      n_obs <- GenEL_get_n_obs(private$.gelptr)
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
      GenEL_omega_hat(private$.gelptr, lambda, t(G), weights)
    },

    #' @description Evaluate the empirical loglikelihood function.
    #'
    #' @template param_G
    #' @template param_weights
    #' @param check_conv If `TRUE`, checks whether the desired `rel_tol` was reached within `max_iter` Newton-Raphson iterations.  If not, returns `-Inf.
    #'
    #' @return The empirical loglikelihood evaluated at `G` (a scalar).
    logel = function(G, weights, check_conv = TRUE) {
      private$check_G(G)
      if (missing(weights) || is.null(weights)) {
        ans <- GenEL_logel(private$.gelptr, t(G), check_conv)
      }
      else {
        ## n_obs <- GenEL_get_n_obs(private$.gelptr)
        ## if (length(weights) != n_obs) {
        ##   stop("Length of `weights` does not equal to the number of obserations.")
        ## }
        ## if (any(weights < 0)) {
        ##   stop("`weights` should contain only non-negative values.")
        ## }
        private$check_weights(weights)
        ans <- GenEL_weighted_logel(private$.gelptr, t(G), weights, check_conv)
      }
      ans
    },

    #' @description Calculate empirical loglikelihood and its derivative with respect to `G`.
    #'
    #' @template param_G
    #' @template param_weights
    #' @param check_conv If `TRUE`, checks whether the desired `rel_tol` was reached within `max_iter` Newton-Raphson iterations.  If not, sets the empirical loglikelihood to `-Inf` and its gradient to a matrix of `NaN`s.
    #'
    #' @return A list with elements `logel` and `grad`, the former being a scalaer and the latter a matrix of size `n_obs x n_eqs`.
    logel_grad = function(G, weights, check_conv = TRUE) {
      private$check_G(G)
      if (missing(weights) || is.null(weights)) {
        ans <- GenEL_logel_grad(private$.gelptr, t(G), check_conv)
      }
      else {
        ## n_obs <- GenEL_get_n_obs(private$.gelptr)
        ## if (length(weights) != n_obs) {
        ##   stop("Length of `weights` does not equal to the number of obserations.")
        ## }
        ## if (any(weights < 0)) {
        ##   stop("`weights` should contain only non-negative values.")
        ## }
        private$check_weights(weights)
        ans <- GenEL_weighted_logel_grad(private$.gelptr, t(G), weights, check_conv)
      }
      ans$grad <- t(ans$grad)
      ans
    }
  )
)
