#include <Rcpp.h>
using namespace Rcpp;
//[[Rcpp::depends("RcppEigen")]]
#include <RcppEigen.h>
using namespace Eigen;
#include "AdjG.h"

// [[Rcpp::export(".adjG")]]
Eigen::MatrixXd adjG(Eigen::MatrixXd G, double a) {
  int nObs = G.cols();
  int nEqs = G.rows();
  Eigen::MatrixXd aG = Eigen::MatrixXd::Zero(nEqs,nObs+1);
  aG.leftCols(nObs) = G;
  adj_G(aG,a);
  return aG;
}
