#ifndef MEANREGLSMODEL_h
#define MEANREGLSMODEL_h

#include <math.h>
#include <Rmath.h>

// subclass: Mean regression location-scale model
class MeanRegLSModel {
private:
    RowVectorXd yXb;   // y'-beta'X
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
    // G = MatrixXd::Zero(nEqs,nObs);
    G = MatrixXd::Zero(nEqs+1,nObs);
    yXb = RowVectorXd::Zero(nObs);
}

// form the G matrix
inline void MeanRegLSModel::evalG(const Ref<const VectorXd>& theta) {
    // extract params for location and scale parts
    VectorXd beta = theta.head(nEqs/2);
    VectorXd gamma = theta.tail(nEqs/2);
    // std::cout << "beta = " << beta << std::endl;
    // std::cout << "gamma = " << gamma << std::endl;
    // location 
    yXb.noalias() = y.transpose() - beta.transpose() * X;
    // std::cout << "yXb = \n" << yXb << std::endl; 
    RowVectorXd gXe = (-2*gamma.transpose() * X).array().exp(); // 1 x nObs
    // std::cout << "gXe = \n" << gXe << std::endl;
    G.block(0,0,nEqs/2,nObs) = X;
    G.block(0,0,nEqs/2,nObs).array().rowwise() *= yXb.array()*gXe.array();
    // std::cout << "G = \n" << G << std::endl;
    // scale 
    RowVectorXd yXb2 = (yXb.array()*yXb.array()).matrix();
    G.block(nEqs/2,0,nEqs/2,nObs) = X;
    G.block(nEqs/2,0,nEqs/2,nObs).array().rowwise() *= yXb2.array()*gXe.array();
    G.block(nEqs/2+1,0,nEqs/2,nObs).array() = yXb2.array()*gXe.array() - 1;
    // std::cout << "y = \n" << y << std::endl;
    std::cout << "G = \n" << G << std::endl;
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
