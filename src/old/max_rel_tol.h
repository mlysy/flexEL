/// @file max_rel_tol.h

#ifndef MAX_REL_TOL_H
#define MAX_REL_TOL_H

#include <RcppEigen.h>

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

}

#endif
