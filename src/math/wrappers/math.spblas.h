
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

template <typename T, template <typename> class Matrix>
struct SpBLAS {
	SpBLAS();
	SpBLAS(Matrix<T> &a);
	~SpBLAS();

	void insert(Matrix<T> &a);
	void insertTransposed(Matrix<T> &a);
	void insert(Matrix_Dense<T> &a, double threshold);

	void extractUpper(Matrix<T> &a);

	// y = alpha * A * x + beta * y
	void apply(Vector_Dense<T> &y, const T &alpha, const T &beta, const Vector_Dense<T> &x);

	// this = A * B
	void multiply(SpBLAS<T, Matrix> &A, SpBLAS<T, Matrix> &B);
	void multiply(SpBLAS<T, Matrix> &A, Matrix_Dense<T> &B);

	// this = A * At
	void AAt(SpBLAS<T, Matrix> &A);

	void transposeTo(SpBLAS<T, Matrix> &A);
	void convertTo(Matrix_Dense<T> &out);

	void solveRowMayor(Matrix_Dense<T> &rhs, Matrix_Dense<T> &solution);
	void solveColMayor(Matrix_Dense<T> &rhs, Matrix_Dense<T> &solution);

	// input = output[start_row:end_row, start_col:end_col]. Start inclusive, end exclusive
	void submatrix(Matrix_Dense<T> &output, esint start_row, esint end_row, esint start_col, esint end_col, bool trans = false, bool conj = false, bool output_force_full = false);
	void submatrix(Matrix_CSR<T>   &output, esint start_row, esint end_row, esint start_col, esint end_col, bool trans = false, bool conj = false, bool output_force_full = false);

	static void submatrix(const Matrix_CSR<T> &input, Matrix_Dense<T> &output, esint start_row, esint end_row, esint start_col, esint end_col, bool trans = false, bool conj = false, bool output_force_full = false);
	static void submatrix(const Matrix_CSR<T> &input, Matrix_CSR<T>   &output, esint start_row, esint end_row, esint start_col, esint end_col, bool trans = false, bool conj = false, bool output_force_full = false);

	Matrix<T> *matrix;

private:
	Matrix_SpBLAS_External_Representation *_spblas;
};

}

#endif /* SRC_MATH_WRAPPERS_MATH_SPBLAS_H_ */
