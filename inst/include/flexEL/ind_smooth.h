/**
 * @file ind_smooth.h
 * 
 * @brief Smoothed indicator function.
 */

#ifndef IND_SMOOTH_H
#define IND_SMOOTH_H

#include <math.h>

namespace flexEL {

  /**
   * @brief Smoothed indicator function.
   * 
   * @param[in] x   The value to be evaluated at. For x > 0, the function returns a number that is 
   *               close to 0; for x < 0, the function returns a number that is close to 1.
   * @param[in] s   A positive number controlling the level of smoothness. The bigger the more smooth.
   */
  inline double ind_smooth(double x, double s) {
    return(1.0/(1.0+exp(x*s)));
  }
  
  /**
   * @brief First derivative of \c ind_smooth.
   * 
   * @param[in] x   The value to be evaluated at. 
   * @param[in] s   A positive number controlling the level of smoothness. The bigger the more smooth.
   */
  inline double ind1_smooth(double x, double s) {
    double sxexp = exp(s*x);
    return(-s*sxexp/((1+sxexp)*(1+sxexp)));
  }

   /**
    * @brief Second derivative of \c ind_smooth.
    * 
    * @param[in] x   The value to be evaluated at. 
    * @param[in] s   A positive number controlling the level of smoothness. The bigger the more smooth.
    */
   inline double ind2_smooth(double x, double s) {
      double sxexp = exp(s*x);
      double sxexpp1 = 1+sxexp;
      return(-s*s*sxexp/(sxexpp1*sxexpp1) + 
             2*s*s*sxexp*sxexp/(sxexpp1*sxexpp1*sxexpp1));
   }

} // namespace flexEL

#endif
