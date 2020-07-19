/**
 * @file dvec_to_barr.h
 * 
 * @brief Convert Eigen VectorXd (double vector) to C++ bool array.
 */

#ifndef DVEC_TO_BARR_H
#define DVEC_TO_BARR_H

// convert Eigen VectorXd dVec to C++ bool array bArr
// bArr must have the same length as dVec
namespace flexEL {

  inline void dVec_to_bArr(Eigen::VectorXd dVec, bool *bArr) {
    // std::cout << "dVec.size() = " << dVec.size() << std::endl;
    for (int ii=0; ii<dVec.size(); ii++) {
      if (dVec(ii) == 0) bArr[ii] = false;
      else bArr[ii] = true;
    }
  }

}

#endif