/// @file utils.h
///
/// Utility functions.  Some of these are provided here rather than as private class methods to simplify testing.

#ifndef FLEXEL_UTILS_H
#define FLEXEL_UTILS_H

#include <Eigen/Dense>

namespace flexEL {
  using namespace Eigen;

  /// Maximum relative error between two vectors.
  ///
  /// Calculates
  /// ```
  /// max |x1 - x2| / (|x1 + x2| + .1)
  /// ```
  /// The constant in the denominator is for situations in which some of thex's are very small in absolute value.
  ///
  /// @param x1[in] First vector.
  /// @param x2[in] Second vector.
  /// @return The maximum relative error.
  inline double max_rel_err(const Ref<const VectorXd>& x1,
                            const Ref<const VectorXd>&x2) {
    return ((x1 - x2).array().abs() /
            ((x1 + x2).array().abs() + 0.1)).maxCoeff();
  }

  /// Smoothed indicator function.
  ///
  /// Calculates `smooth_ind = 1/(1 + exp(s * (eps1 - eps2)))`, where `eps1` is a scalar and `eps2` is a vector.
  ///
  /// @param smooth_ind[out] Vector in which to store the result.
  /// @param eps1[in] First residual.
  /// @param eps2[in] Second residual.
  /// @param s[in] Smoothing parameter.
  inline void smooth_indicator(Ref<VectorXd> smooth_ind,
			       double eps1,
			       const Ref<const VectorXd>& eps2,
			       double s) {
    smooth_ind = (1.0 + (s * (eps1 - eps2.array())).exp()).inverse();
    return;
  }

  /// Calculate the adjusted G matrix.
  ///
  /// Returns the adjusted `G` matrix given a `G` matrix of dimension `n_eqs x n_obs+1`.  That is, it modifies the last column to be the artificial observation rather than creating a new matrix.
  ///
  /// Reference: J. Chen, A. M. Variyath, and B. Abraham. "Adjusted empirical likelihood and its properties".  Journal of Computational and Graphical Statistics, 17(2):426–443, 2008.
  ///
  /// @param[in/out] G A numeric matrix of size `n_eqs x (n_obs + 1)`.  Only the last column is overwritten by the calculation.
  /// @param[in] a Scalar tuning parameter for the adjustment.
  inline void adj_G(Ref<MatrixXd> G, double a) {
    // std::cout << "adj_G: G = \n" << G << std::endl;
    int n_obs = G.cols()-1;
    // Eigen::VectorXd gbar = 1.0/n_obs*G.rowwise().sum();
    // std::cout << "gbar = " << gbar.transpose() << std::endl;
    // G.rightCols(1) = -a*gbar;
    G.rightCols(1) = -(a/n_obs) * G.leftCols(n_obs).rowwise().sum();
    return;
  }

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

}

#endif
