library(testthat)

if(!identical(Sys.getenv("NOT_CRAN"), "true")) {
  set.seed(230)
}

test_check("flexEL")
