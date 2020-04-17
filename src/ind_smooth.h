/**
 * @file IndSmooth.h
 * 
 * @brief Smoothed indicator function.
 */

#ifndef INDSMOOTH_h
#define INDSMOOTH_h

#include <math.h>

namespace flexEL {

  inline double ind_smooth(double x, double s) {
    return(1.0/(1.0+exp(x*s)));
  }
  
  inline double ind1_smooth(double x, double s) {
    double sxexp = exp(s*x);
    return(-s*sxexp/((1+sxexp)*(1+sxexp)));
  }

} // namespace flexEL

// #include <RcppEigen.h>
// using namespace Eigen;

//// [[Rcpp::depends(RcppEigen)]]

// inline VectorXd ind_smooth(VectorXd x, VectorXd s) {
//   return(1.0/(1.0+(x.array()*s.array()).exp()));
// }
// 
// inline VectorXd ind1_smooth(VectorXd x, VectorXd s) {
//   ArrayXd sxexp = (s.array()*x.array()).exp();
//   return(-s.array()*sxexp/((1+sxexp)*(1+sxexp)));
// }

#endif