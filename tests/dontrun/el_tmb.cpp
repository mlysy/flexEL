/// @file el_tmb.cpp
///
/// @brief TMB wrapper for `flexEL::GenEL<Type>.logel()`.
///
/**
    While `flexEL` is fully templated, we wouldn't want to use autodiff on the Newton-Raphson solver directly.  Rather, we would use `flexEL::GenEL<Type>.logel_grad()`.  Here's how this would work:
    
    1. `TMB_ATOMIC_VECTOR_FUNCTION` will create a function `logel_atomic()` with a custom reverse-mode rule for us.

    - The forward pass (`ATOMIC_DOUBLE`) has input and output of type `CppAD::vector<double>`:

        - The input (`tx`) should be a flattened version of `G`, prepended with the row/col dimensions and any other information, e.g., support correction, number of NR iterations, etc.

	- The output (`ty`) would normally be the result of `flexEL::GenEL<double>.logel(G)`, which is a scalar.  However, to efficiently compute the reverse-mode algorithm we actually need the intermediate results obtainted by calling `flexEL::GenEL<double>.logel_full()`, which returns `logel` as scalar and `omega` and `lambda` as intermediate vectors.  These vectors then need to be concatenated into a vector.  Thus, the `OUTPUT_DIM` macro needs to be set to the sum of the elements in `G`, `omega`, and `lambda`.

    - The backward pass (`ATOMIC_REVERSE`) consists of the body of a function with inputs and outputs all being of type `CppAD::vector<Type>`.

        -  The inputs are `tx`, `ty`, and `py`.  The first two correspond to the same quantities as in `ATOMIC_DOUBLE`, but having been converted to type `CppAD::vector<Type>`.  `py` is a vector of the same length as `ty`.

        - The output is `px`, a vector of the same length as `tx`.

	- The body of the backward pass first unpacks `G` from `tx`, and `omega` and `lambda` from `ty`.  Then it uses `flexEL::GenEL<Type>.logel_grad()` to compute `dldG`, the gradient of `logel_atomic()` with respect to `G`.  Finally, we set 
            ```
            px[i+n_nondiff] = d/dG[i] sum(logel(G) * py[0]) =  
                  = dldG[i] * py[0],
            ```
	    where `n_nondiff` is the number of nondifferentiable arguments prepended to the flattened version of `G` (i.e., dimensions of `G` and NR control parameters).

    2.  `logel_atomic()` is a function which takes vector inputs and outputs, each of which contain the flattened version of a matrix (`G` or `dldG`) along with some elements (corresponding to the nondifferentiable inputs) which will never get used.  So the final `logel()` function is just a wrapper which does the arguent packing and unpacking for us.

*/

#include "flexEL/gen_el.h"
#include <TMB.hpp>
#include <Eigen/Dense>
#include <Eigen/Core>
using namespace Eigen;

// typedefs
template <class Type>
using cRefMatrix_t = const Ref <const Matrix<Type, Dynamic, Dynamic> >;
template <class Type>
using MatrixX_t = Matrix<Type, Dynamic, Dynamic>;
template <class Type>
using VectorX_t = Matrix<Type, Dynamic, 1>;

/// Copy contiguous elements of a `CppAD::vector<Type>` to a `VectorX_t<Type>`.
///
/// The intented usage is something like
/// ```
/// const VectorX_t<Type> y = copy_vector<Type>(x)
/// ```
/// And this _does_ make a copy of the data: see <https://stackoverflow.com/questions/52261389/how-to-convert-an-stdvector-to-a-matrix-in-eigen>.
///
/// @param[in] x A `CppAD::vector<Type>` object.
/// @param[in] start The starting index.
/// @param[in] size The number of elements to map.
/// @return A `VectorX_t<Type>` consisting of the elements `x[start], ..., x[start + size-1]`.
template <class Type>
const VectorX_t<Type> copy_vector(const CppAD::vector<Type>& x,
			    int start, int size) {
  return Map<const VectorX_t<Type>,Unaligned>(x.data()+start, size);
}

/// Copy contiguous elements of a `CppAD::vector<Type>` to a `MatrixX_t<Type>`.
///
/// @param[in] x A `CppAD::vector<Type>` object.
/// @param[in] start The starting index.
/// @param[in] rows The number of rows in the output matrix.
/// @param[in] cols The number of columns in the output matrix.
/// @return A `Matrix_t<Type>` object with dimensions `rows x cols` and consisting of the elements `x[start], ..., x[start + (row * cols)-1]`.
template <class Type>
const MatrixX_t<Type> copy_matrix(const CppAD::vector<Type>& x,
			    int start, int rows, int cols) {
  return Map<const MatrixX_t<Type>,Unaligned>(x.data()+start, rows, cols);
}


template <class Type>
const MatrixX_t<Type> tx_to_G(const CppAD::vector<Type>& tx) {
  int n_row = CppAD::Integer(tx[0]);
  int n_col = CppAD::Integer(tx[1]);
  return copy_matrix<Type>(tx, 2, n_row, n_col);
}

/// Vectorized wrapper to `flexEL::GenEL<Type>.logel_full()`.
///
/// @param[out] ty A vector containing the result, `omega`, and `lambda`.
/// @param[in] tx A vector containing `G.rows()`, `G.cols()`, and the flattened `G`.
template <class Type>
void logel_vector(CppAD::vector<Type>& ty, const CppAD::vector<Type>& tx) {
  // unpack from tx
  int n_eqs = int(tx[0]);
  int n_obs = int(tx[1]);
  const MatrixX_t<Type> G = copy_matrix(tx, 2, n_eqs, n_obs);
  // weights input needed for logel_full()
  VectorX_t<Type> weights = VectorX_t<Type>::Ones(n_obs);
  // outputs
  VectorX_t<Type> lambda(n_eqs);
  VectorX_t<Type> omega(n_obs);
  // calculations
  flexEL::GenEL<Type> gel(n_obs, n_eqs);
  ty[0] = gel.logel_full(omega, lambda, G, weights);
  // pack into ty
  // NOTE: this makes a copy of the Eigen objects,
  // see <https://stackoverflow.com/questions/26094379/typecasting-eigenvectorxd-to-stdvector>
  VectorX_t<Type>::Map(&ty[1], omega.size()) = omega;
  VectorX_t<Type>::Map(&ty[1+omega.size()], lambda.size()) = lambda;
  return;
}

TMB_ATOMIC_VECTOR_FUNCTION(
			   // ATOMIC_NAME
			   logel_atomic
			   ,
			   // OUTPUT_DIM
			   CppAD::Integer(tx[0] * tx[1]),
			   // ATOMIC_DOUBLE
			   logel_vector<double>(ty, tx),
			   // ATOMIC_REVERSE
			   // unpack inputs from tx
			   MatrixX_t<Type> G = tx_to_G(tx);
			   int n_eqs = G.rows();
			   int n_obs = G.cols();
			   // unpack inputs from ty
			   int start = 1;
			   const VectorX_t<Type> omega = copy_vector<Type>(ty, start, n_obs);
			   start += n_obs;
			   const VectorX_t<Type> lambda = copy_vector<Type>(ty, start, n_eqs);
			   // calculate gradient
			   MatrixX_t<Type> dldG(n_eqs, n_obs);
			   flexEL::GenEL<Type> gel(n_obs, n_eqs);
			   gel.logel_grad(dldG, omega, lambda, Type(n_obs));
			   // compute the reverse-mode rule
			   px[0] = Type(0.);
			   px[1] = Type(0.);
			   for(int jj=0; jj<n_obs; jj++) {
			     for(int ii=0; ii<n_eqs; ii++) {
			       px[jj*n_eqs + ii + 2] = dldG(ii,jj) * py[0];
			       // px[jj*n_eqs + ii + 2] = Type(2) * py[0];
			     }
			   }
			   // for(int ii=2; ii<tx.size(); ii++) {
			   //   px[ii] = py[0];
			   // }
)

/// Wrapper to logel_atomic() for Eigen types.
///
/// @param G Matrix of size `n_eqs x n_obs`.
/// @return The return value of `flexEL::GenEL<Type>.logel()`.
template<class Type>
Type logel(cRefMatrix_t<Type>& G) {
  CppAD::vector<Type> tx(G.size() + 2);
  int n_row = G.rows();
  int n_col = G.cols();
  tx[0] = Type(n_row);
  tx[1] = Type(n_col);
  for(int jj=0; jj<n_col; jj++) {
    for(int ii=0; ii<n_row; ii++) {
      tx[jj*n_row + ii + 2] = G(ii,jj);
    }
  }
  return logel_atomic(tx)[0];
}

template<class Type>
Type objective_function<Type>::operator() ()
{
  PARAMETER_MATRIX(G);
  Type f = logel<Type>(G);
  return f;
}
