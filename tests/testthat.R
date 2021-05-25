library(testthat)

if(!identical(Sys.getenv("NOT_CRAN"), "true")) {
  set.seed(123)
}

test_check("flexEL")
