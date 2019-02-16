/**
 * @file SortOrder.h
 * 
 * @brief Find the order of elements if a vector is sorted (ascendingly).
 */

#ifndef SORTORDER_h
#define SORTORDER_h

#include <RcppEigen.h>
using namespace Eigen;
// [[Rcpp::depends(RcppEigen)]]

// Modifed from: https://stackoverflow.com/a/12399290/5186269

// 
/**
 * @file       SortOrder.h
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

// return indices after sorting 
// TODO: this could be a void function and pass the inds as an argument too 
inline VectorXi sort_inds(const VectorXd &vec) {
    int n = vec.size();
    VectorXi inds(n);
    for (int ii=0; ii<n; ii++) {
        inds(ii) = ii;
    }
    std::sort(inds.data(), inds.data()+n, compare_acc_vec <decltype(vec)> (vec));
    return inds;
}

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