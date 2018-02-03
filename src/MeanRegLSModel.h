#ifndef MEANREGLSMODEL_h
#define MEANREGLSMODEL_h

#include <math.h>
#include <Rmath.h>

// subclass: Mean regression location-scale model
class MeanRegLSModel {
private:
    RowVectorXd yXb;   // y'-beta'X
    RowVectorXd yXblg; // log(yXb)-gamma'X
    MatrixXd tG;
protected:
    VectorXd y;
    MatrixXd X;
    int nObs, nEqs;
    MatrixXd G;
public:
    MeanRegLSModel(const Ref<const VectorXd>& _y, const Ref<const MatrixXd>& _X,
                   void* params);
    void evalG(const Ref<const VectorXd>& theta);
};

// constructor
inline MeanRegLSModel::MeanRegLSModel(const Ref<const VectorXd>& _y,
                                      const Ref<const MatrixXd>& _X,
                                      void* params) {
    y = _y;
    X = _X;
    nObs = y.size();
    nEqs = X.rows() * 2; // nObs x 2 for the mean and location function
    G = MatrixXd::Zero(nEqs,nObs);
    tG = MatrixXd::Zero(nObs, nEqs);
    yXb = RowVectorXd::Zero(nObs);
    yXblg = RowVectorXd::Zero(nObs);
}

// form the G matrix
inline void MeanRegLSModel::evalG(const Ref<const VectorXd>& theta) {
    // extract params for location and scale parts
    VectorXd beta = theta.head(nEqs/2);
    VectorXd gamma = theta.tail(nEqs/2);
    // location (mean) part
    yXb.noalias() = y.transpose() - beta.transpose() * X;
    tG.block(0,0,nObs,nEqs/2) = X.transpose();
    tG.block(0,0,nObs,nEqs/2).array().colwise() *= yXb.transpose().array();
    G.block(0,0,nObs,nEqs/2) = tG.block(0,0,nObs,nEqs/2).transpose();
    // scale (variance) part
    yXblg.noalias() = yXb.array().abs().log().matrix() - gamma.transpose() * X;
    tG.block(0,nEqs,nObs,nEqs/2) = X.transpose();
    tG.block(0,nEqs,nObs,nEqs/2).array().colwise() *= yXblg.transpose().array();
}

#endif
