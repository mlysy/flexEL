#ifndef QUANTREGMODEL_h
#define QUANTREGMODEL_h

#include <math.h>
#include <Rmath.h>

// subclass: Mean regression
class QuantRegModel {
private:
    double rho_alpha(double u, double alpha); // TODO: not needed?
    double phi_alpha(double u, double alpha); 
protected:
    VectorXd y;
    MatrixXd X;
    double alpha; 
    int nObs, nEqs;
    MatrixXd G;
public:
    QuantRegModel(const Ref<const VectorXd>& _y, const Ref<const MatrixXd>& _X,
                  void* params); // constructor non-censor
    void evalG(const Ref<const VectorXd>& beta);
    MatrixXd getG(); // TODO: should prob move to EL since duplicate for MR and QR
};

// constructor
inline QuantRegModel::QuantRegModel(const Ref<const VectorXd>& _y, 
                                    const Ref<const MatrixXd>& _X,
                                    void* params) {
    y = _y;
    X = _X; 
    alpha = *(double*)(params); 
    nObs = y.size();
    nEqs = X.rows(); // X gets passed as p x nObs matrix
    G = MatrixXd::Zero(nEqs,nObs);
}

// revised L1 loss function for quantile regression
inline double QuantRegModel::rho_alpha(double u, double alpha) {
    return(u * (alpha - (u <= 0)));
}

// 1st derivative of rho_alpha
inline double QuantRegModel::phi_alpha(double u, double alpha) {
    return((u <= 0) - alpha);
}

// form the G matrix
inline void QuantRegModel::evalG(const Ref<const VectorXd>& beta) {
    for(int ii=0; ii<y.size(); ii++) {
        this->G.col(ii) = phi_alpha(y(ii)-X.col(ii).transpose()*beta, alpha)*X.col(ii);
    }
}

// get function for G matrix
inline MatrixXd QuantRegModel::getG() {
    return(G);
}

#endif