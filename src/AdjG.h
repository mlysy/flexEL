#ifndef ADJG_h
#define ADJG_h

#include <RcppEigen.h>

// returns the adjusted G matrix given a G matrix of dimension nEqs x nObs
inline Eigen::MatrixXd adj_G(Eigen::MatrixXd G, double a) {
  int nObs = G.cols();
  int nEqs = G.rows();
  Eigen::VectorXd gbar = 1.0/nObs*G.rowwise().sum();
  Eigen::VectorXd gadd = -a*gbar;
  Eigen::MatrixXd aG(nEqs,nObs+1);
  aG.block(0,0,nEqs,nObs) = G;
  aG.block(0,nObs,nEqs,1) = gadd;
  return(aG);
}

#endif