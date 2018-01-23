

class QuantRegModel {
 private:
  double alpha;
  VectorXd y;
  MatrixXd X;
  int nObs, nEq;
 public:
  MatrixXd G;
  QuantRegModel(VectorXd _y, MatrixXd _X, void *par);
  ~QuantRegModel;
  void evalG(VectorXd theta);
};

inline QuantRegModel::QuantRegModel(VectorXd _y, MatrixXd _X, void *par) {
  y = _y; // make sure it's a COPY, not just pointer allocation
  X = _X;
  nObs = y.length(); // or something
  nEq = X.row();
  G = MatrixXd::Zero(nObs, nEq); // or however you allocate Eigen memory
  alpha = &((double*)par)
}

// theta is beta.
inline void QuantRegModel::evalG(VectorXd theta) {
}

inline QuantRegModel::~QuantRegModel {}

// ---------------------------------------------------------------------------------------------------------------------------------------------------------

template <class elMod>
class InnerEL : elMod {
  // make sure G and evalG are correctly inherited
  
};
