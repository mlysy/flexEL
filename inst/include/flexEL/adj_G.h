/**
 * @file adj_G.h
 * 
 * @brief Calculate the adjusted G matrix.
 */

#ifndef ADJ_G_H
#define ADJ_G_H

#include <RcppEigen.h>

// returns the adjusted G matrix given a G matrix of dimension nEqs x nObs+1
// Note: it modifies the last row to be the artificial obs rather than adding one row
namespace flexEL {

  /**
   * @brief Calculate the adjusted G matrix given a G matrix.
   * 
   * J. Chen, A. M. Variyath, and B. Abraham. Adjusted empirical likelihood and its properties.
   *    Journal of Computational and Graphical Statistics, 17(2):426–443, 2008.
   * 
   * @param G   A numeric matrix.
   * @param a   Tuning parameter for the adjustment.
   */
  inline void adj_G(Eigen::MatrixXd &G, double a) {
    // std::cout << "adj_G: G = \n" << G << std::endl;
    int nObs = G.cols()-1;
    Eigen::VectorXd gbar = 1.0/nObs*G.rowwise().sum();
    // std::cout << "gbar = " << gbar.transpose() << std::endl;
    G.rightCols(1) = -a*gbar;
  }

}

#endif
