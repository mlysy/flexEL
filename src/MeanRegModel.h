#ifndef MEANREGMODEL_h
#define MEANREGMODEL_h

#include <math.h>
#include <Rmath.h>

// subclass: Mean regression
class MeanRegModel {
private:
    RowVectorXd yXb;
    // MatrixXd tG;
protected:
    VectorXd y;
    MatrixXd X;
    int nObs, nEqs;
    MatrixXd G;
public:
    MeanRegModel(const Ref<const VectorXd>& _y, const Ref<const MatrixXd>& _X,
                 void* params);
    void evalG(const Ref<const VectorXd>& beta);
};

// constructor
inline MeanRegModel::MeanRegModel(const Ref<const VectorXd>& _y,
                                  const Ref<const MatrixXd>& _X,
                                  void* params) {
    y = _y;
    X = _X;
    nObs = y.size();
    nEqs = X.rows(); // X gets passed as p x nObs matrix
    G = MatrixXd::Zero(nEqs,nObs);
    // tG = MatrixXd::Zero(nObs, nEqs);
    yXb = RowVectorXd::Zero(nObs);
}

// form the G matrix
inline void MeanRegModel::evalG(const Ref<const VectorXd>& beta) {
    // yXb.noalias() = y.transpose() - beta.transpose() * X;
    // tG = X.transpose();
    // tG.array().colwise() *= yXb.transpose().array();
    // G = tG.transpose();
    yXb.noalias() = y.transpose() - beta.transpose() * X;
    G = X;
    G.array().rowwise() *= yXb.array();
    // std::cout << "G = " << G << std::endl;
}

#endif
