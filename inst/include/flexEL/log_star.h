/// @file log_star.h

#ifndef LOG_STAR_H
#define LOG_STAR_H

namespace flexEL {

  /// Modified log function for the Newton-Raphson dual problem.
  ///
  /// The modified log (or log-star) function is defined as `log(x)` for `x > trunc`, and `ax^2 + bx + c` for `x < trunc`.
  ///
  /// @param[in] x Value at which to evaluate the function.
  /// @param[in] trunc Truncation value.
  /// @param[in] a Polynomial coefficient of degree two.
  /// @param[in] b Polynomial coefficient of degree one.
  /// @param[in] c Polynomial coefficient of degree zero.
  ///
  /// @return Value of the log-star function.
  inline double log_star(double x,
			 double trunc, double a, double b, double c) {
    if(x >= trunc) {
      return(log(x));
    } else {
      return((a*x + b)*x + c);
    }
  }

  /// First derivative of the log-star function.
  ///
  /// @param[in] x Value at which to evaluate the function.
  /// @param[in] trunc Truncation value.
  /// @param[in] a Polynomial coefficient of degree two.
  /// @param[in] b Polynomial coefficient of degree one.
  ///
  /// @return Value of the first derivative of the log-star function.
  inline double log_star1(double x,
			  double trunc, double a, double b) {
    if(x >= trunc) {
      return(1.0/x);
    } 
    else {
      return(2.0*a*x + b);
    }
  }

  /// Second derivative of the log-star function.
  ///
  /// @param[in] x Value at which to evaluate the function.
  /// @param[in] trunc Truncation value.
  /// @param[in] a Polynomial coefficient of degree two.
  ///
  /// @return Value of the first derivative of the log-star function.
  inline double log_star2(double x, double trunc, double a) {
    if(x >= trunc) {
      return(-1.0/(x*x));
    } else {
      return(2.0*a);
    }
  }
  

} // end namespace flexEL

#endif
