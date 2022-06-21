#--- API examples for flexEL ---------------------------------------------------

#' Function to calculate `G` for EL regression.
#'
mrG <- function(beta, y, X) {
  # code goes here
}


#--- Basic API -----------------------------------------------------------------

el <- GenEL$new(n_obs, n_eq)
mr_el <- function(beta) {
  el$logel(mrG(beta, y, X))
}

# can do something similar for G and dG/dbeta.

#--- Slightly less basic, many for multiple usage ------------------------------

MREL <- R6Class(
  classname = "MREL",
  inherit = "GenEL",
  public = list(
    initialize = function(X, y)
  )
  )

#--- desired API for mean regression -------------------------------------------

el <- MREL$new(y, X, Z = NULL, delta = NULL, el_opts = NULL)

to_theta = function(beta, gamma = NULL) {}

to_coef = function(theta) {
  # returns either beta or a list with elements beta and gamma
}
