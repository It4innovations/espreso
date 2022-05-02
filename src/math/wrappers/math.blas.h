
#ifndef SRC_MATH_WRAPPERS_MATH_BLAS_H_
#define SRC_MATH_WRAPPERS_MATH_BLAS_H_

#include "math/primitives/vector_dense.h"
#include "math/primitives/vector_sparse.h"
#include "math/primitives/matrix_dense.h"
#include "math/primitives/matrix_csr.h"
#include "math/primitives/matrix_ijv.h"

namespace espreso {
namespace math {

	// x = y
	template <typename T>
	void copy(const esint size, T *x, const int incX, const T *y, const int incY);

	// x *= alpha
	template <typename T>
	void scale(const esint size, const T &alpha, T *x, const int incX);

	// x += alpha * y
	template <typename T>
	void add(const esint size, T *x, const int incX, const T &alpha, const T *y, const int incY);

	template <typename T>
	T dot(const esint size, const T *x, const int incX, const T *y, const int incY);

	template <typename T>
	T norm(const esint size, const T *x, const int incX);

	// y = alpha * A * x + beta * y
	template <typename T>
	void apply(Vector_Dense<T> &y, const T &alpha, const Matrix_Dense<T> &a, const T &beta, const Vector_Dense<T> &x);

	// y = alpha * At * x + beta * y
	template <typename T>
	void applyT(Vector_Dense<T> &y, const T &alpha, const Matrix_Dense<T> &a, const T &beta, const Vector_Dense<T> &x);
}
}


namespace espreso {
namespace math {

	template <typename T> void copy(Vector_Dense<T>  &x, const Vector_Dense<T>  &y) { copy(x.size           , x.vals, 1, y.vals, 1); }
	template <typename T> void copy(Vector_Sparse<T> &x, const Vector_Sparse<T> &y) { copy(x.nnz            , x.vals, 1, y.vals, 1); }
	template <typename T> void copy(Matrix_Dense<T>  &x, const Matrix_Dense<T>  &y) { copy(x.nrows * x.ncols, x.vals, 1, y.vals, 1); }
	template <typename T> void copy(Matrix_CSR<T>    &x, const Matrix_CSR<T>    &y) { copy(x.nnz            , x.vals, 1, y.vals, 1); }
	template <typename T> void copy(Matrix_IJV<T>    &x, const Matrix_IJV<T>    &y) { copy(x.nnz            , x.vals, 1, y.vals, 1); }

	template <typename T> void scale(const T &alpha, Vector_Dense<T>  &x) { scale(x.size           , alpha, x.vals, 1); }
	template <typename T> void scale(const T &alpha, Vector_Sparse<T> &x) { scale(x.nnz            , alpha, x.vals, 1); }
	template <typename T> void scale(const T &alpha, Matrix_Dense<T>  &x) { scale(x.nrows * x.ncols, alpha, x.vals, 1); }
	template <typename T> void scale(const T &alpha, Matrix_CSR<T>    &x) { scale(x.nnz            , alpha, x.vals, 1); }
	template <typename T> void scale(const T &alpha, Matrix_IJV<T>    &x) { scale(x.nnz            , alpha, x.vals, 1); }

	template <typename T> void add(Vector_Dense<T>  &x, const T &alpha, const Vector_Dense<T>  &y) { add(x.size           , x.vals, 1, alpha, y.vals, 1); }
	template <typename T> void add(Vector_Sparse<T> &x, const T &alpha, const Vector_Sparse<T> &y) { add(x.nnz            , x.vals, 1, alpha, y.vals, 1); }
	template <typename T> void add(Matrix_Dense<T>  &x, const T &alpha, const Matrix_Dense<T>  &y) { add(x.nrows * x.ncols, x.vals, 1, alpha, y.vals, 1); }
	template <typename T> void add(Matrix_CSR<T>    &x, const T &alpha, const Matrix_CSR<T>    &y) { add(x.nnz            , x.vals, 1, alpha, y.vals, 1); }
	template <typename T> void add(Matrix_IJV<T>    &x, const T &alpha, const Matrix_IJV<T>    &y) { add(x.nnz            , x.vals, 1, alpha, y.vals, 1); }

	template <typename T> T dot(const Vector_Dense<T>  &x, const Vector_Dense<T>  &y) { return dot(x.size, x.vals, 1, y.vals, 1); }
	template <typename T> T dot(const Vector_Sparse<T> &x, const Vector_Sparse<T> &y) { return dot(x.nnz , x.vals, 1, y.vals, 1); }

	template <typename T> T norm(const Vector_Dense<T>  &x) { return norm(x.size, x.vals, 1); }
	template <typename T> T norm(const Vector_Sparse<T> &x) { return norm(x.nnz , x.vals, 1); }
}
}

#endif /* SRC_MATH_WRAPPERS_MATH_BLAS_H_ */
