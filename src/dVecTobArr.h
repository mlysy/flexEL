#ifndef DVECTOBARR_h
#define DVECTOBARR_h

// convert Eigen VectorXd dVec to C++ bool array bArr
// bArr must have the same length as dVec
inline void dVec_to_bArr(Eigen::VectorXd dVec, bool *bArr) {
  // std::cout << "dVec.size() = " << dVec.size() << std::endl;
  for (int ii=0; ii<dVec.size(); ii++) {
    if (dVec(ii) == 0) bArr[ii] = false;
    else bArr[ii] = true;
  }
}

#endif