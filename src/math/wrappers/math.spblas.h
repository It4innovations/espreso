
#ifndef SRC_MATH_WRAPPERS_MATH_SPBLAS_H_
#define SRC_MATH_WRAPPERS_MATH_SPBLAS_H_

#include "math/primitives/vector_dense.h"
#include "math/primitives/vector_sparse.h"
#include "math/primitives/matrix_dense.h"
#include "math/primitives/matrix_csr.h"
#include "math/primitives/matrix_ijv.h"

namespace espreso {
namespace math {

	// utility functions allowing the Intel inspector-executor model
	template <typename T> void commit(Matrix_Dense<T> &x);
	template <typename T> void commit(Matrix_CSR<T> &x);
	template <typename T> void commit(Matrix_IJV<T> &x);

	template <typename T> void free(Matrix_Dense<T> &x);
	template <typename T> void free(Matrix_CSR<T> &x);
	template <typename T> void free(Matrix_IJV<T> &x);

	// y = alpha * A * x + beta * y
	template <typename T> void apply(Vector_Dense<T> &y, const T &alpha, const Matrix_CSR<T>   &a, const T &beta, const Vector_Dense<T> &x);
	template <typename T> void apply(Vector_Dense<T> &y, const T &alpha, const Matrix_IJV<T>   &a, const T &beta, const Vector_Dense<T> &x);

	// input = output[start_row:end_row, start_col:end_col]. Start inclusive, end exclusive
	template <typename T> void submatrix(const Matrix_CSR<T> &input, Matrix_Dense<T> &output, esint start_row, esint end_row, esint start_col, esint end_col, bool trans = false, bool conj = false, bool output_force_full = false);
	template <typename T> void submatrix(const Matrix_CSR<T> &input, Matrix_CSR<T>   &output, esint start_row, esint end_row, esint start_col, esint end_col, bool trans = false, bool conj = false, bool output_force_full = false);

}
}

#endif /* SRC_MATH_WRAPPERS_MATH_SPBLAS_H_ */
