/// @file block_outer.h

#ifndef BLOCK_OUTER_H
#define BLOCK_OUTER_H

#include <RcppEigen.h>

namespace flexEL {

  using namespace Eigen;

  /// For an `m x N` matrix `G = [g1 | ... | gN]`, calculates the `m x (m*N)` matrix `GGt = [g1 g1' | ... | gN gN']`.
  ///
  /// @param[out] GGt Matrix of size `m x (m*N)` in which the calculated result will be saved.
  /// @param[in] G Matrix of size `m x N`.
  inline void block_outer(Ref<MatrixXd> GGt, const Ref<const MatrixXd>& G) {
    int n_eqs = G.rows();
    int n_obs = G.cols();
    // for each row of G, compute outer product and store as block and put into GGt
    for(int ii=0; ii<n_obs; ii++) {
      GGt.block(0,ii*n_eqs,n_eqs,n_eqs).noalias() = G.col(ii) * G.col(ii).transpose();
    }
    return;
  }

}

#endif
