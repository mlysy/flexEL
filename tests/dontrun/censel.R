#--- Prototype for R CensEL class ----------------------------------------------

require(flexEL)
require(R6)

# wrapper to GenEL that allows us to access its internal GEL.

GenEL_ex <- R6::R6Class(
  classname = "GenEL_ex",
  inherit = GenEL,
  active = list(
    gel = function() {
      super$get_gel()
    }
  )
)

gel <- GenEL_ex$new(10, 3)
