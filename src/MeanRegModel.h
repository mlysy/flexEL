#ifndef MEANREGMODEL_h
#define MEANREGMODEL_h

// #include <math.h>
// #include <Rmath.h>

class MeanRegModel {
private:
    RowVectorXd yXb; 
    MatrixXd tG;
protected:
    VectorXd y;
    MatrixXd X;
    int nObs, nEqs;
    MatrixXd G;
public:
    // MeanRegModel(const Ref<const VectorXd>& _y, const Ref<const MatrixXd>& _X,
    //              void* params); // old ctor 
    // MeanRegModel(); // default ctor -- it shouldn't have one actually 
    void setData(const Ref<const VectorXd>& _y, const Ref<const MatrixXd>& _X,
                 void* params); // set data with default ctor
    // for location linear regression models
    void evalG(const Ref<const VectorXd>& beta); 
    void setG(const Ref<const MatrixXd>& _G); 
    MatrixXd getG();
};

/*
// constructor: old
inline MeanRegModel::MeanRegModel(const Ref<const VectorXd>& _y,
                                  const Ref<const MatrixXd>& _X,
                                  void* params) {
    y = _y;
    X = _X;
    nObs = y.size();
    nEqs = X.rows(); // X gets passed as p x nObs matrix
    G = MatrixXd::Zero(nEqs,nObs);
    tG = MatrixXd::Zero(nObs, nEqs);
    yXb = RowVectorXd::Zero(nObs);
}
*/

// setData (with default ctor)
inline void MeanRegModel::setData(const Ref<const VectorXd>& _y,
                                  const Ref<const MatrixXd>& _X,
                                  void* params) {
    y = _y;
    X = _X;
    nObs = y.size();
    nEqs = X.rows(); // X gets passed as p x nObs matrix
    G = MatrixXd::Zero(nEqs,nObs);
    tG = MatrixXd::Zero(nObs, nEqs);
    yXb = RowVectorXd::Zero(nObs);
}

// form the G matrix for location linear regression model 
inline void MeanRegModel::evalG(const Ref<const VectorXd>& beta) {
    yXb.noalias() = y.transpose() - beta.transpose() * X;
    tG = X.transpose();
    tG.array().colwise() *= yXb.transpose().array();
    G = tG.transpose();
    // Note: rowwise is slower in general 
    // yXb.noalias() = y.transpose() - beta.transpose() * X;
    // G = X;
    // G.array().rowwise() *= yXb.array();
    // std::cout << "G = " << G << std::endl;
}

// set function for G matrix
inline void MeanRegModel::setG(const Ref<const MatrixXd>& _G) {
    G = _G; 
}

// get function for G matrix
inline MatrixXd MeanRegModel::getG() {
    return(G);
}

#endif
