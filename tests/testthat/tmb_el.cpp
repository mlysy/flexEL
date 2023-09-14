/// @file tmb_el.cpp
///
/// @brief Test TMB wrapper for `flexEL::GenEL<Type>.logel()`.

#include <TMB.hpp>
#include "flexEL/tmb_el.hpp"

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_INTEGER(max_iter);
  DATA_SCALAR(rel_tol);
  DATA_INTEGER(supp_adj);
  PARAMETER_MATRIX(G);
  Type f = flexEL::logel<Type>(G, max_iter, rel_tol, supp_adj);
  return f;
}
