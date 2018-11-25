/**
 * @file MeanRegModel.h
 * 
 * @brief A base class for template classes InnerEL and InnerELC.
 */

#ifndef MEANREGMODEL_h
#define MEANREGMODEL_h

// #include <math.h>
// #include <Rmath.h>

/**
 * @file       MeanRegModel.h
 * 
 * @class      MeanRegModel
 * 
 * @brief      A base class for template classes InnerEL and InnerELC to calculate estimating equations for mean regressions.
 */
class MeanRegModel {
private:

  RowVectorXd yXb;
  RowVectorXd eZg;
  RowVectorXd yXbeZg;
  RowVectorXd yXbeZg2;
  MatrixXd tG;

protected:

  int nObs, nEqs, nBet, nGam;
  VectorXd y;
  MatrixXd X;
  MatrixXd Z;
  MatrixXd G;

public:

  // constructors
  /**
  * @brief Default constructor for MeanRegModel.
  */
  MeanRegModel();

  /**
   * @brief Constructor for MeanRegModel with dimensions as inputs.
   * 
   * @param nObs    Number of observations.
   * @param nEqs    Number of estimating equations.
   */
  MeanRegModel(int nObs, int nEqs);

  // set data
  /**
   * @brief Set data for mean regression location model.
   * 
   * @param y      Responses of length \code{nObs}.
   * @param X      Covariate matrix of dimension \code{nBet} x \code{nObs}.
   */
  void setData(const Ref<const VectorXd>& y,
               const Ref<const MatrixXd>& X);

  /**
   * @brief Set data for mean regression location-scale model.
   * 
   * @param y      Responses of length \code{nObs}.
   * @param X      Covariate matrix of dimension \code{nBet} x \code{nObs}.
   * @param Z      Covariate matrix of dimension \code{nGam} x \code{nObs}.
   */
  void setData(const Ref<const VectorXd>& y,
               const Ref<const MatrixXd>& X,
               const Ref<const MatrixXd>& Z);

  // evalutate G matrix
  /**
   * @brief Evaluate G matrix for mean regression location model.
   * 
   * @param beta     Coefficient vector in linear location function.
   */
  void evalG(const Ref<const VectorXd>& beta);

  /**
   * @brief Evaluate G matrix for mean regression location-scale model.
   * 
   * @param beta     Coefficient vector of length \code{nBet} in linear location function.
   * @param gamma    Coefficient vector of length \code{nGam} in exponential scale function.
   * @param sig2     Scale parameter in scale function.
   */
  void evalG(const Ref<const VectorXd>& beta,
             const Ref<const VectorXd>& gamma,
             const double& sig2);
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
  y = _y;
  X = _X;
  nObs = y.size();
  nBet = X.rows(); // X gets passed as nBet x nObs matrix
  nGam = 0;
  nEqs = nBet;
  G = MatrixXd::Zero(nEqs,nObs);
  tG = MatrixXd::Zero(nObs, nEqs);
  yXb = RowVectorXd::Zero(nObs);
}

// setData location-scale model (with default ctor)
inline void MeanRegModel::setData(const Ref<const VectorXd>& _y,
                                  const Ref<const MatrixXd>& _X,
                                  const Ref<const MatrixXd>& _Z) {
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
  yXb.noalias() = y.transpose() - beta.transpose() * X;
  tG = X.transpose();
  tG.array().colwise() *= yXb.transpose().array();
  G = tG.transpose();
}

// form the G matrix for location-scale linear regression model
inline void MeanRegModel::evalG(const Ref<const VectorXd>& beta,
                                const Ref<const VectorXd>& gamma,
                                const double &sig2) {
  eZg.array() = (-gamma.transpose()*Z).array().exp();
  yXbeZg.array() = (y.transpose()-beta.transpose()*X).array() * eZg.array();
  yXbeZg2.array() = yXbeZg.array()*yXbeZg.array();
  tG.block(0,0,nObs,nBet) = X.transpose();
  tG.block(0,0,nObs,nBet).array().colwise() *= yXbeZg.transpose().array() * eZg.transpose().array();
  tG.block(0,nBet,nObs,nGam) = Z.transpose();
  tG.block(0,nBet,nObs,nGam).array().colwise() *= (1.0-yXbeZg2.transpose().array());
  tG.rightCols(1).array() = 1/sig2*yXbeZg2.transpose().array()-1;
  G = tG.transpose();
}

#endif
