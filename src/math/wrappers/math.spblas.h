
#ifndef SRC_MATH_WRAPPERS_MATH_SPBLAS_H_
#define SRC_MATH_WRAPPERS_MATH_SPBLAS_H_

#include "math/primitives/vector_dense.h"
#include "math/primitives/vector_sparse.h"
#include "math/primitives/matrix_dense.h"
#include "math/primitives/matrix_csr.h"
#include "math/primitives/matrix_csc.h"
#include "math/primitives/matrix_ijv.h"

namespace espreso {

struct Matrix_SpBLAS_External_Representation;

template <template <typename, typename, typename> class Matrix, typename T, typename I = int>
struct SpBLAS {
	using MatrixType = Matrix<T, I, cpu_allocator>;

	SpBLAS();
	SpBLAS(MatrixType &a);
	~SpBLAS();

	void insert(MatrixType &a);

	// y = alpha * A * x + beta * y
	void apply(Vector_Dense<T, I> &y, const T &alpha, const T &beta, const Vector_Dense<T, I> &x);

	// input = output[start_row:end_row, start_col:end_col]. Start inclusive, end exclusive
	void submatrix(Matrix_Dense<T, I> &output, I start_row, I end_row, I start_col, I end_col, bool trans = false, bool conj = false, bool output_force_full = false);
	void submatrix(Matrix_CSR<T, I>   &output, I start_row, I end_row, I start_col, I end_col, bool trans = false, bool conj = false, bool output_force_full = false);

	static void submatrix(const Matrix_CSR<T, I> &input, Matrix_Dense<T, I> &output, I start_row, I end_row, I start_col, I end_col, bool trans = false, bool conj = false, bool output_force_full = false);
	static void submatrix(const Matrix_CSR<T, I> &input, Matrix_CSR<T, I>   &output, I start_row, I end_row, I start_col, I end_col, bool trans = false, bool conj = false, bool output_force_full = false);

	MatrixType *matrix;

private:
	Matrix_SpBLAS_External_Representation *_spblas;
};

}

#include "math.spblas.hpp"

#endif /* SRC_MATH_WRAPPERS_MATH_SPBLAS_H_ */
