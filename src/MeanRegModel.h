#ifndef MEANREGMODEL_h
#define MEANREGMODEL_h

// #include <math.h>
// #include <Rmath.h>

class MeanRegModel {
private:
  
  RowVectorXd yXb;
  RowVectorXd eZg;
  RowVectorXd yXbeZg;
  RowVectorXd yXbeZg2; 
  MatrixXd tG;
  
protected:
  
  VectorXd y;
  MatrixXd X;
  MatrixXd Z;
  int nObs, nEqs, nBet, nGam, nQts; // TODO: nQts is not used in mean reg
  MatrixXd G;
  
public:
  
  // constructors
  /**
  * @brief Default constructor for MeanRegModel.
  */
  MeanRegModel();
  
  /**
   * @brief Constructor for MeanRegModel with dimensions as inputs.
   * @param _nObs    Number of observations.
   * @param _nEqs    Number of estimating equations.
   */
  MeanRegModel(int _nObs, int _nEqs);
  
  // set data
  /**
   * @brief Set data for mean regression location model.
   * @param _y      Responses of length nObs.
   * @param _X      Covariate matrix of dimension \code{nBet} x \code{nObs}.
   */
  void setData(const Ref<const VectorXd>& _y, 
               const Ref<const MatrixXd>& _X);
               // void* params);
  
  /**
   * @brief Set data for mean regression location-scale model.
   * @param _y      Responses of length nObs.
   * @param _X      Covariate matrix of dimension \code{nBet} x \code{nObs}.
   * @param _Z      Covariate matrix of dimension \code{nGam} x \code{nObs}.
   */
  void setData(const Ref<const VectorXd>& _y, 
               const Ref<const MatrixXd>& _X,
               const Ref<const MatrixXd>& _Z);
               // void* params);
            
  // evalutate G matrix   
  /**
   * @brief Evaluate G matrix for mean regression location model.
   * @param beta     Coefficient vector in linear location function.
   */
  void evalG(const Ref<const VectorXd>& beta); 
  
  /**
   * @brief Evaluate G matrix for mean regression location-scale model.
   * @param beta     Coefficient vector of length \code{nBet} in linear location function.
   * @param gamma    Coefficient vector of length \code{nGam} in exponential scale function.
   * @param sig2     Scale parameter in scale function.
   */
  void evalG(const Ref<const VectorXd>& beta, 
             const Ref<const VectorXd>& gamma, 
             const double& sig2,
             const Ref<const VectorXd>& dummy);
  
  // setter and getter for G matrix
  void setG(const Ref<const MatrixXd>& _G); 
  MatrixXd getG();
};

// default ctor 
inline MeanRegModel::MeanRegModel(){}

// ctor
inline MeanRegModel::MeanRegModel(int _nObs, int _nEqs) {
  nObs = _nObs;
  nEqs = _nEqs; // X gets passed as nBet x nObs matrix
  G = MatrixXd::Zero(nEqs,nObs);
}

// setData location model (with default ctor)
inline void MeanRegModel::setData(const Ref<const VectorXd>& _y,
                                  const Ref<const MatrixXd>& _X) {
                                  // void* params) {
  y = _y;
  X = _X;
  nObs = y.size();
  nBet = X.rows(); // X gets passed as nBet x nObs matrix
  nGam = 0; 
  nEqs = nBet; 
  G = MatrixXd::Zero(nBet,nObs);
  tG = MatrixXd::Zero(nObs, nBet);
  yXb = RowVectorXd::Zero(nObs);
}

// setData location-scale model (with default ctor)
inline void MeanRegModel::setData(const Ref<const VectorXd>& _y,
                                  const Ref<const MatrixXd>& _X,
                                  const Ref<const MatrixXd>& _Z) {
                                  // void* params) {
  y = _y;
  X = _X;
  Z = _Z;
  nObs = y.size();
  nBet = X.rows(); // X gets passed as nBet x nObs matrix
  nGam = Z.rows(); // Z gets passed as nGam x nObs matrix
  nEqs = nBet+nGam+1; 
  G = MatrixXd::Zero(nEqs,nObs);
  tG = MatrixXd::Zero(nObs, nEqs);
  eZg = RowVectorXd::Zero(nObs);
  yXbeZg = RowVectorXd::Zero(nObs);
  yXbeZg2 = RowVectorXd::Zero(nObs);
}

// form the G matrix for location linear regression model 
inline void MeanRegModel::evalG(const Ref<const VectorXd>& beta) {
  // std::cout << "evalG: beta = " << beta.transpose() << std::endl;
  yXb.noalias() = y.transpose() - beta.transpose() * X;
  tG = X.transpose();
  tG.array().colwise() *= yXb.transpose().array();
  G = tG.transpose();
}

// form the G matrix for location-scale linear regression model 
inline void MeanRegModel::evalG(const Ref<const VectorXd>& beta, 
                                const Ref<const VectorXd>& gamma,
                                const double &sig2,
                                const Ref<const VectorXd>& dummy) {
  eZg.array() = (-gamma.transpose()*Z).array().exp();
  yXbeZg.array() = (y.transpose()-beta.transpose()*X).array() * eZg.array();
  yXbeZg2.array() = yXbeZg.array()*yXbeZg.array();
  tG.block(0,0,nObs,nBet) = X.transpose();
  tG.block(0,0,nObs,nBet).array().colwise() *= yXbeZg.transpose().array() * eZg.transpose().array();
  tG.block(0,nBet,nObs,nGam) = Z.transpose();
  // tG.block(0,nBet,nObs,nGam).array().colwise() *= yXbeZg2.transpose().array();
  tG.block(0,nBet,nObs,nGam).array().colwise() *= (1.0-yXbeZg2.transpose().array());
  tG.rightCols(1).array() = 1/sig2*yXbeZg2.transpose().array()-1;
  G = tG.transpose();
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
