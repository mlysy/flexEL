#' R6 class for EL regression where response is under right censoring.
#' 
#' A general EL object.
#' @export
CensEL <- R6::R6Class(
  
  classname = "CensEL",
  
  private = list(
    .CEL = NULL
  ),
  
  active = list(),
  
  public = list(
    
    #' @description Create a new CensEL object.
    #' @param n_obs Number of observations.
    #' @param n_eqs Number of (moment constraint) equations.
    #' @return A `CensEL` object.
    initialize = function(n_obs, n_eqs) {
      private$.CEL <- GenEL_ctor(n_obs, n_eqs)
      private$.lambda0 <- rep(0, n_eqs)
    }
    
  )
  
)