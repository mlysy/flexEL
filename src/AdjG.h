#ifndef ADJG_h
#define ADJG_h

#include <RcppEigen.h>

// returns the adjusted G matrix given a G matrix of dimension nEqs x nObs+1
// Note: it modifies the last row to be the artificial obs rather than adding one row
inline void adj_G(Eigen::MatrixXd &G, double a) {
  int nObs = G.cols()-1;
  // int nEqs = G.rows();
  Eigen::VectorXd gbar = 1.0/nObs*G.rowwise().sum();
  // std::cout << "gbar = " << gbar.transpose() << std::endl;
  // Eigen::VectorXd gadd = -a*gbar;
  G.rightCols(1) = -a*gbar;
  // Eigen::MatrixXd aG(nEqs,nObs+1);
  // aG.block(0,0,nEqs,nObs) = G;
  // aG.block(0,nObs,nEqs,1) = gadd;
  // return(aG);
}
// TODO, just modify the last row of G matrix instead of adding one row

#endif