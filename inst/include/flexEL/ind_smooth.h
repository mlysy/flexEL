/**
 * @file ind_smooth.h
 * 
 * @brief Smoothed indicator function.
 */

#ifndef IND_SMOOTH_H
#define IND_SMOOTH_H

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

#endif
