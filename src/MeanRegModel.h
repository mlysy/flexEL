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

  RowVectorXd yXb_; // y - X*beta
  RowVectorXd eZg_; // exp(-Z*gamma)
  RowVectorXd yXbeZg_; // (y - X*beta)*exp(-Z*gamma)
  RowVectorXd yXbeZg2_; // (y - X*beta)^2*exp(-2*Z*gamma)
  MatrixXd tG_; // t(G)
  int nBet_, nGam_; // dimensions of beta and gamma
  VectorXd y_; // vector of responses
  MatrixXd X_; // covariate matrix used in location function
  MatrixXd Z_; // covariate matrix used in scale function

protected:

  int nObs_; /**< number of observations (number of columns of G) */
  int nEqs_; /**< number of estimating equations (number of rows of G) */
  MatrixXd G_; /**< matrix of estimating equations of dimension <code>nEqs</code> x <code>nObs</code>*/

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
   * @param y      Responses of length <code>nObs</code>.
   * @param X      Covariate matrix of dimension <code>nBet</code> x <code>nObs</code>.
   */
  void setData(const Ref<const VectorXd>& y,
               const Ref<const MatrixXd>& X);

  /**
   * @brief Set data for mean regression location-scale model.
   * 
   * @param y      Responses of length <code>nObs</code>.
   * @param X      Covariate matrix of dimension <code>nBet</code> x <code>nObs</code>.
   * @param Z      Covariate matrix of dimension <code>nGam</code> x <code>nObs</code>.
   */
  void setData(const Ref<const VectorXd>& y,
               const Ref<const MatrixXd>& X,
               const Ref<const MatrixXd>& Z);

  // evalutate G matrix
  /**
   * @brief Evaluate G matrix for mean regression location model.
   * 
   * @param beta     Coefficient vector of length <code>nObs</code> in linear location function.
   */
  void evalG(const Ref<const VectorXd>& beta);

  /**
   * @brief Evaluate G matrix for mean regression location-scale model.
   * 
   * @param beta     Coefficient vector of length <code>nBet</code> in linear location function.
   * @param gamma    Coefficient vector of length <code>nGam</code> in exponential scale function.
   * @param sig2     Scale parameter in scale function.
   */
  void evalG(const Ref<const VectorXd>& beta,
             const Ref<const VectorXd>& gamma,
             const double& sig2);
};

// default ctor
inline MeanRegModel::MeanRegModel(){}

// ctor
inline MeanRegModel::MeanRegModel(int nObs, int nEqs) {
  nObs_ = nObs;
  nEqs_ = nEqs; // X gets passed as nBet x nObs matrix
  // G_ = MatrixXd::Zero(nEqs_,nObs_);
}

// setData location model (with default ctor)
inline void MeanRegModel::setData(const Ref<const VectorXd>& y,
                                  const Ref<const MatrixXd>& X) {
  y_ = y;
  X_ = X;
  nObs_ = y.size();
  nBet_ = X.rows(); // X gets passed as nBet x nObs matrix
  nGam_ = 0;
  nEqs_ = nBet_;
  // G_ = MatrixXd::Zero(nEqs_,nObs_);
  // tG_ = MatrixXd::Zero(nObs_, nEqs_);
  // yXb_ = RowVectorXd::Zero(nObs_);
}

// setData location-scale model (with default ctor)
inline void MeanRegModel::setData(const Ref<const VectorXd>& y,
                                  const Ref<const MatrixXd>& X,
                                  const Ref<const MatrixXd>& Z) {
  y_ = y;
  X_ = X;
  Z_ = Z;
  nObs_ = y.size();
  nBet_ = X.rows(); // X gets passed as nBet x nObs matrix
  nGam_ = Z.rows(); // Z gets passed as nGam x nObs matrix
  nEqs_ = nBet_+nGam_+1;
  // G_ = MatrixXd::Zero(nEqs_,nObs_);
  // tG_ = MatrixXd::Zero(nObs_, nEqs_);
  // eZg_ = RowVectorXd::Zero(nObs_);
  // yXbeZg_ = RowVectorXd::Zero(nObs_);
  // yXbeZg2_ = RowVectorXd::Zero(nObs_);
}

// form the G matrix for location linear regression model
inline void MeanRegModel::evalG(const Ref<const VectorXd>& beta) {
  // yXb_.noalias() = y_.transpose() - beta.transpose() * X_;
  yXb_ = y_.transpose() - beta.transpose() * X_;
  tG_ = X_.transpose();
  tG_.array().colwise() *= yXb_.transpose().array();
  G_ = tG_.transpose();
}

// form the G matrix for location-scale linear regression model
inline void MeanRegModel::evalG(const Ref<const VectorXd>& beta,
                                const Ref<const VectorXd>& gamma,
                                const double &sig2) {
  eZg_.array() = (-gamma.transpose()*Z_).array().exp();
  yXbeZg_.array() = (y_.transpose()-beta.transpose()*X_).array() * eZg_.array();
  yXbeZg2_.array() = yXbeZg_.array()*yXbeZg_.array();
  
  tG_ = MatrixXd::Zero(nObs_, nEqs_); // NEW: DEC 25
  tG_.block(0,0,nObs_,nBet_) = X_.transpose();
  tG_.block(0,0,nObs_,nBet_).array().colwise() *= yXbeZg_.transpose().array() * eZg_.transpose().array();
  tG_.block(0,nBet_,nObs_,nGam_) = Z_.transpose();
  tG_.block(0,nBet_,nObs_,nGam_).array().colwise() *= (1.0-yXbeZg2_.transpose().array());
  tG_.rightCols(1).array() = 1/sig2*yXbeZg2_.transpose().array()-1;
  G_ = tG_.transpose();
}

#endif
