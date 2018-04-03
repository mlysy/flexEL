#ifndef MEANREGLSMODEL_h
#define MEANREGLSMODEL_h

#include <math.h>
#include <Rmath.h>

// subclass: Mean regression location-scale model
class MeanRegLSModel {
private:
    RowVectorXd yXb;   // y'-beta'X
    RowVectorXd gZe; // exp(-2 * gamma'Z)
    VectorXd WW; // yXb * gZe
    MatrixXd tG;
    int nBeta;
    int nGamma;
    VectorXd y;
    MatrixXd tX;
    MatrixXd X;
    MatrixXd tZ;
    MatrixXd Z;
    VectorXd beta;
    VectorXd gamma;
protected:
    int nObs, nEqs;
    MatrixXd G;
public:
    MeanRegLSModel(const Ref<const VectorXd>& _y,
		   const Ref<const MatrixXd>& _X,
		   const Ref<const MatrixXd>& _Z,
                   void* params);
    void evalG(const Ref<const VectorXd>& theta);
};

// constructor
inline MeanRegLSModel::MeanRegLSModel(const Ref<const VectorXd>& _y,
                                      const Ref<const MatrixXd>& _X,
				      const Ref<const MatrixXd>& _Z,
                                      void* params) {
    y = _y;
    X = _X; // on regular scale
    tX = X.transpose();
    Z = _Z;
    tZ = Z.transpose();
    nBeta = tX.rows();
    nGamma = tZ.rows();
    nObs = y.size();
    nEqs = nBeta + nGamma + 1; // nObs x 2 for the mean and location function
    // G = MatrixXd::Zero(nEqs,nObs);
    G = MatrixXd::Zero(nEqs,nObs);
    tG = G.transpose();
    yXb = RowVectorXd::Zero(nObs);
    gZe = RowVectorXd::Zero(nObs);
    WW = VectorXd::Zero(nObs);
    beta = VectorXd::Zero(nBeta);
    gamma = VectorXd::Zero(nGamma);
}

// form the G matrix
inline void MeanRegLSModel::evalG(const Ref<const VectorXd>& theta) {
    // extract params for location and scale parts
    beta = theta.head(nBeta);
    gamma = theta.tail(nGamma);
    // std::cout << "beta = " << beta << std::endl;
    // std::cout << "gamma = " << gamma << std::endl;
    // location
    yXb.noalias() = - beta.transpose() * tX;
    yXb.array() += y.transpose().array();
    //yXb.noalias() = (y.array() - (beta.transpose() * tX).array()).transpose();
    // std::cout << "yXb = \n" << yXb << std::endl;
    gZe.noalias() = -2.0 * gamma.transpose() * tZ;
    gZe = gZe.array().exp();
    // std::cout << "gXe = \n" << gXe << std::endl;
    tG.leftCols(nBeta) = X;
    WW = yXb.array()*gZe.array();
    tG.leftCols(nBeta).array().colwise() *= WW.array();
    // std::cout << "G = \n" << G << std::endl;
    // scale
    WW.array() *= yXb.array();
    tG.block(0, nBeta, nObs, nGamma) = Z;
    tG.block(0, nBeta, nObs, nGamma).array().colwise() *= WW.array();
    tG.col(nEqs-1) = WW.array() - 1.0;
    G = tG.transpose();
    // std::cout << "y = \n" << y << std::endl;
    // std::cout << "G = \n" << G << std::endl;
}

// inline void MeanRegLSModel::evalG(const Ref<const VectorXd>& theta) {
//     // extract params for location and scale parts
//     VectorXd beta = theta.head(nEqs/2);
//     VectorXd gamma = theta.tail(nEqs/2);
//     std::cout << "beta = " << beta << std::endl;
//     std::cout << "gamma = " << gamma << std::endl;
//     // location 
//     yXb.noalias() = y.transpose() - beta.transpose() * X;
//     RowVectorXd gX = (gamma.transpose() * X).array(); // 1 x nObs
//     std::cout << "gX = \n" << gX << std::endl;
//     RowVectorXd egX = gX.array().exp().matrix(); // 1 x nObs
//     std::cout << "egX = \n" << egX << std::endl; 
//     std::cout << "yXb = \n" << yXb << std::endl; 
//     RowVectorXd yXbegX = (yXb.array() * egX.array()).matrix();  // 1 x nObs
//     std::cout << "bXbegX = \n" << yXbegX << std::endl;
//     double c = yXbegX.dot(yXbegX) / (double)nObs - 1.0; // TODO: needed ??
//     std::cout << "c = " << c << std::endl;
//     RowVectorXd gX2e = (-2*gX.array()).exp().matrix(); // 1 x nObs
//     G.block(0,0,nEqs/2,nObs) = X;
//     G.block(0,0,nEqs/2,nObs).array().colwise() *= ((2*c*gX2e.array() + 1) * yXb.array()).transpose();
//     // scale 
//     G.block(nEqs/2,0,nEqs/2,nObs) = X;
//     G.block(nEqs/2,0,nEqs/2,nObs).array().colwise() *= 
//         (c*(yXb.array()*yXb.array())*gX2e.array()).transpose();
//     std::cout << "G = \n" << G << std::endl;
// }

// inline void MeanRegLSModel::evalG(const Ref<const VectorXd>& theta) {
//     // extract params for location and scale parts
//     VectorXd beta = theta.head(nEqs/2);
//     VectorXd gamma = theta.tail(nEqs/2);
//     // location (mean) 
//     yXb.noalias() = y.transpose() - beta.transpose() * X;
//     G.block(0,0,nEqs/2,nObs) = X;
//     ArrayXd g = yXb.array() * ((-2*gamma.transpose() * X).array().exp());
//     G.block(0,0,nEqs/2,nObs).array().colwise() *= g;
//     // scale (variance) 
//     // yXblg.noalias() = yXb.array().abs().log().matrix() - gamma.transpose() * X;
//     // tG.block(0,nEqs/2,nObs,nEqs/2) = X.transpose();
//     // tG.block(0,nEqs/2,nObs,nEqs/2).array().colwise() *= yXblg.transpose().array();
//     // G.block(nEqs/2,0,nEqs/2,nObs) = tG.block(0,nEqs/2,nObs,nEqs/2).transpose();
//     g = yXb.transpose().array() * g; 
//     double h = g.sum()/nObs - 1;
//     std::cout << h << std::endl;
//     G.block(nEqs/2,0,nEqs/2,nObs) = h*X;
//     G.block(nEqs/2,0,nEqs/2,nObs).array().colwise() *= g;
// }

#endif
