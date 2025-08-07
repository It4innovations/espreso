
#ifndef SRC_MATH_WRAPPERS_MATH_LAPACK_H_
#define SRC_MATH_WRAPPERS_MATH_LAPACK_H_

#include "math/primitives/matrix_dense.h"
#include "math/primitives/vector_dense.h"

namespace espreso {
namespace math {
namespace lapack {

template <typename T, typename I>
void solve_sym_upper(Matrix_Dense<T, I> &A, Matrix_Dense<T, I> &rhs);

template <typename T, typename I>
void solve_general(Matrix_Dense<T, I> &A, Matrix_Dense<T, I> &rhs);

template <typename T, typename I>
void get_eig_sym(Matrix_Dense<T, I> &A, Vector_Dense<T, I> &values);

template <typename T, typename I>
void get_eig_sym(Matrix_Dense<T, I> &A, Vector_Dense<T, I> &values, Matrix_Dense<T, I> &vectors);

template <typename T, typename I>
void get_eig_sym(Matrix_Dense<T, I> &A, Vector_Dense<T, I> &values, I begin, I end);

template <typename T, typename I>
void get_eig_sym(Matrix_Dense<T, I> &A, Vector_Dense<T, I> &values, Matrix_Dense<T, I> &vectors, I begin, I end);

template <typename T, typename I>
void get_svd(Matrix_Dense<T, I> &A, Vector_Dense<double> &s, Matrix_Dense<double> &U, Matrix_Dense<double> &V);

// input = output[start_row:end_row, start_col:end_col]. Start inclusive, end exclusive
template <typename T, typename I>
void submatrix(const Matrix_Dense<T, I> &input, Matrix_Dense<T, I> &output, I start_row, I end_row, I start_col, I end_col);

}
}
}

#endif /* SRC_MATH_WRAPPERS_MATH_LAPACK_H_ */
