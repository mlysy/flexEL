/**
 * @file CensELExports.cpp
 * 
 * @brief Rcpp wrappers for GenEL class.
 */

#include <Rcpp.h>
#include <RcppEigen.h>
#include "inner_elc_new.h"

//[[Rcpp::depends("RcppEigen")]]

/// Construct a CensEL object.
///
/// Instantiates a `CensEL` object on the C++ side, wraps it in an `Rcpp::XPtr`, and returns the corresponding `externalptr` on the R side.
///
/// @param[in] n_obs    Number of observations.
/// @param[in] n_eqs    Number of estimating equations.
/// @return An `externalptr` pointing to the CensEL object.
///
// [[Rcpp::export]]
SEXP CensEL_ctor(int n_obs, int n_eqs) {
  flexEL::CensEL *CEL = new flexEL::CensEL(n_obs, n_eqs);
  Rcpp::XPtr<flexEL::CensEL> pGEL(CEL, true);
  return pGEL;
}




