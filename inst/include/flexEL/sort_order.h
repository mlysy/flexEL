
// Modifed from: https://stackoverflow.com/a/12399290/5186269

#ifndef SORT_ORDER_H
#define SORT_ORDER_H

#include <RcppEigen.h>

using namespace Eigen;

// [[Rcpp::depends(RcppEigen)]]

namespace flexEL {

  /**
   * @file       sort_order.h
   *
   * @class      compare_acc_vec
   *
   * @brief      A helper class for obtaining indices corresponding to ascendingly sorted values in a vector
   */
  template <typename T>
  class compare_acc_vec {
      const T& vec;
  public:
      compare_acc_vec(const T& vec): vec(vec) { }
      bool operator () (size_t lind, size_t rind) const {
          return vec[lind] < vec[rind]; // < to sort ascendingly; > to sort descendingly
      }
  };

  // TODO: this could be a void function and pass the inds as an argument too 
  /**
   * @brief Find the indices of elements in a vector if it is sorted ascendingly.
   * 
   * @param vec A numeric vector.
   * 
   * @return A vector of integers (the indicies of the elements in \c vec).
   */
  inline VectorXi sort_inds(const VectorXd &vec) {
      int n = vec.size();
      VectorXi inds(n);
      for (int ii=0; ii<n; ii++) {
          inds(ii) = ii;
      }
      std::sort(inds.data(), inds.data()+n, compare_acc_vec <decltype(vec)> (vec));
      return inds;
  }

} // namespace flexEL

// Original version: 
// Note: have to change return type to vector<int> if want to convert to
//       Eigen Vector types, i.e., VectorXi
// template <typename T>
// vector<size_t> sort_inds(const vector<T> &vec) {
//     vector<size_t> inds(vec.size());
//     iota(inds.begin(), inds.end(), 0);
//     sort(inds.begin(), inds.end(), compare_acc_vec <decltype(vec)> (vec));
//     return inds;
// }

#endif
