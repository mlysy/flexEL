

GenEL <- R6::R6Class(
  
  classname = "GenEL",
  
  private = list(
    .y = NULL,
    .X = NULL,
    
    .max_iter = 100,
    .rel_tol = 1e-7,
    .supp_corr = TRUE,
    .lambda0 = 0
  ),
  
  active = list(
    
    .y = function(value) {
      if (missing(value)) private$.y
      else private$.y <- value
    },
    
    .X = function(value) {
      if (missing(value)) private$.X
      else private$.X <- value
    }
  ),
  
  public = list(
    
    initialize = function(y, X) {
      private$.y <- y
      private$.X <- X
    },
    
    set_opts = function(max_iter, rel_tol, supp_corr, lambda0) {
      .max_iter <- max_iter
      .rel_tol <- rel_tol
      .supp_corr <- supp_corr
      .lambda0 <- lambda0
    },
    
    logel = function(G, ...) {
      logEL(G, ...)
    },
    
    lambda_nr = function(G, ...) {
      lambdaNR(G, ...)
    }
  )
)