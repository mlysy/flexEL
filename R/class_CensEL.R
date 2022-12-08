#' R6 Class for censored EL models.
#'
#' @export
CensEL <- R6::R6Class(

                classname = "CensEL",

                private = list(
                  .celptr = NULL,
                  .max_iter = NULL,
                  .abs_tol = NULL,
                  .smooth_s = NULL,

                  #' @description Check that `delta` is well specified.
                  #' @param delta Vector of censoring indicators.
                  check_delta = function(delta) {
                    if(length(delta) != self$n_obs) {
                      stop("`delta` must be of length `n_obs`.")
                    }
                    if(anyNA(is.logical(delta))) {
                      stop("`delta` cannot be coerced to a logical vector.")
                    }
                  },

                  #' @description Check that `epsilon` is well specified.
                  #' @param epsilon Vector of residuals.
                  check_epsilon = function(epsilon) {
                    if(length(epsilon) != self$n_obs) {
                      stop("`epsilon` must be of length `n_obs`.")
                    }
                    if(!is.numeric(epsilon)) {
                      stop("`epsilon` must be a numeric vector.")
                    }
                  },

                  #' @description Check that `omega` is well specified.
                  #' @param omega Vector of residuals.
                  check_omega = function(omega) {
                    if(!(length(omega) %in% (self$n_obs + 0:1))) {
                      stop("`omega` must be of length `n_obs + supp_adj`.")
                    }
                    if(!is.numeric(omega) || !all(omega > 0)) {
                      stop("`omega` must be a numeric vector.")
                    }
                  }


                ),

                active = list(
                  #' @field n_obs Number of observations.
                  n_obs = function() {
                    CensEL_get_n_obs(private$.celptr)
                  },

                  #' @field n_eqs Number of equations.
                  n_eqs = function() {
                    CensEL_get_n_eqs(private$.celptr)
                  },

                  #' @field max_iter Maximum number of EM iterations (positive integer).
                  max_iter = function(value) {
                    if (missing(value)) {
                      private$.max_iter
                    } else if (!is.numeric(value) | value <= 0) {
                      stop("`max_iter` must be a positive integer.")
                    } else {
                      private$.max_iter <- value
                      CensEL_set_max_iter(private$.celptr, as.integer(private$.max_iter))
                    }
                  },

                  #' @field abs_tol Absolute tolerance for terminating the EM algorithm (positive number).
                  abs_tol = function(value) {
                    if (missing(value)) {
                      private$.abs_tol
                    } else if (!is.numeric(value) | value <= 0) {
                      stop("`abs_tol` must be a positive number.")
                    } else {
                      private$.abs_tol <- value
                      CensEL_set_abs_tol(private$.celptr, private$.abs_tol)
                    }
                  },

                  #' @field smooth_s Indicator smoothing parameter.  A positive number, where `smooth_s = Inf` gives the regular indicator function.
                  smooth_s = function(value) {
                    if (missing(value)) {
                      private$.smooth_s
                    } else if (!is.numeric(value) | value <= 0) {
                      stop("`smooth_s` must be a positive number.")
                    } else {
                      private$.smooth_s <- value
                      CensEL_set_smooth(private$.celptr, private$.smooth_s)
                    }
                  }

                ),

                public = list(

                  #' @description Create a new `CensEL` object.
                  #' @template param_gel
                  #' @template param_smooth_s
                  #' @template param_max_iter_em
                  #' @template param_abs_tol
                  initialize = function(gel,
                                        smooth_s = 10, max_iter = 10,
                                        abs_tol = 1e-7) {
                    n_obs <- gel$n_obs
                    n_eqs <- gel$n_eqs
                    private$.celptr <- CensEL_ctor(n_obs = n_obs,
                                                   n_eqs = n_eqs)
                    self$smooth_s <- smooth_s
                    self$max_iter <- max_iter
                    self$abs_tol <- abs_tol
                  },

                  #' @description Set multiple options together.
                  #'
                  #' @template param_smooth_s
                  #' @template param_max_iter_em
                  #' @template param_abs_tol
                  #'
                  #' @details Any number of these arguments can be supplied.  Those that are missing are not modified.
                  set_opts = function(smooth_s, max_iter, abs_tol) {
                    if(!missing(smooth_s)) {
                      self$smooth_s <- smooth_s
                    }
                    if(!missing(max_iter)) {
                      self$max_iter <- max_iter
                    }
                    if(!missing(abs_tol)) {
                      self$abs_tol <- abs_tol
                    }
                  },

                  #' @description Calculate the expected weights in the EM algorithm.
                  #'
                  #' @template param_delta
                  #' @template param_epsilon
                  #' @template param_omega
                  #' @template param_gel
                  #'
                  #' @return Weight vector of length `n_obs + supp_adj`.
                  #'
                  #' @details When `supp_adj = TRUE` it is assumed that `delta[n_obs+1] = FALSE` and `epsilon[n_obs+1] = -Inf`.
                  expected_weights = function(delta, epsilon, omega) {
                    private$check_delta(delta)
                    private$check_epsilon(epsilon)
                    private$check_omega(omega)
                    delta <- as.logical(delta)
                    omega <- omega/sum(omega)
                    if(length(omega) == (self$n_obs + 1)) {
                      delta <- c(delta, FALSE)
                      epsilon <- c(epsilon, -Inf)
                    }
                    CensEL_expected_weights(private$.celptr,
                                            delta, epsilon, omega)
                  }
                )
              )
