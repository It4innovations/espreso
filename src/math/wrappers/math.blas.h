
#ifndef SRC_MATH_WRAPPERS_MATH_BLAS_H_
#define SRC_MATH_WRAPPERS_MATH_BLAS_H_

#include "math/primitives/vector_dense.h"
#include "math/primitives/vector_sparse.h"
#include "math/primitives/matrix_dense.h"
#include "math/primitives/matrix_csr.h"
#include "math/primitives/matrix_ijv.h"

#include <complex>

namespace espreso {
namespace math {
namespace blas {

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

	template <typename T>
	T norm(const esint size, const std::complex<T> *x, const int incX);

	// y = alpha * A * x + beta * y
	template <typename T>
	void apply(Vector_Dense<T> &y, const T &alpha, const Matrix_Dense<T> &a, const T &beta, const Vector_Dense<T> &x);

	// y = alpha * At * x + beta * y
	template <typename T>
	void applyT(Vector_Dense<T> &y, const T &alpha, const Matrix_Dense<T> &a, const T &beta, const Vector_Dense<T> &x);

	template <typename T>
	void AAt(const Matrix_Dense<T> &A, Matrix_Dense<T> &AAt);

	// C = alpha * op(A) * op(B) + beta C
	template <typename T>
	void multiply(T alpha, const Matrix_Dense<T> &A, const Matrix_Dense<T> &B, T beta, Matrix_Dense<T> &C, bool transA = false, bool transB = false);
}
}
}

#endif /* SRC_MATH_WRAPPERS_MATH_BLAS_H_ */
