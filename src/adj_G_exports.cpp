/**
 * @file adj_G_exports.cpp
 * 
 * @brief Export the function to calculate adjusted G matrix.
 */

#include <Rcpp.h>
using namespace Rcpp;
//[[Rcpp::depends("RcppEigen")]]
#include <RcppEigen.h>
using namespace Eigen;
#include "adj_G.h"

// [[Rcpp::export(".adjG")]]
Eigen::MatrixXd adjG(Eigen::MatrixXd G, double a) {
  int nObs = G.cols();
  int nEqs = G.rows();
  Eigen::MatrixXd aG = Eigen::MatrixXd::Zero(nEqs,nObs+1);
  aG.leftCols(nObs) = G;
  flexEL::adj_G(aG,a);
  return aG;
}
