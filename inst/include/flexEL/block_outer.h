/**
 * @file block_outer.h
 * 
 * @brief For an (m x N) matrix G = [g1 ... gN], calculates the (m x mN) matrix GGt = [g1 g1' ... gN gN'].
 */

#ifndef BLOCK_OUTER_H
#define BLOCK_OUTER_H

#include <RcppEigen.h>

// for an (m x N) matrix G = [g1 ... gN], calculates the (m x mN) matrix
// GGt = [g1 g1' ... gN gN']
namespace flexEL {

  inline void block_outer(Eigen::MatrixXd &GGt, const Eigen::Ref<const Eigen::MatrixXd>& G) {
    int nEqs = G.rows();
    int nObs = G.cols();
    // for each row of G, compute outer product and store as block and put into GGt
    for(int ii=0; ii<nObs; ii++) {
      GGt.block(0,ii*nEqs,nEqs,nEqs).noalias() = G.col(ii) * G.col(ii).transpose();
    }
    return;
  }

}

#endif