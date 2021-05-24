/// @file adj_G.h

#ifndef ADJ_G_H
#define ADJ_G_H

#include <RcppEigen.h>

namespace flexEL {
  using namespace Eigen;

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

}

#endif
