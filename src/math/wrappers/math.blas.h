
#ifndef SRC_MATH_WRAPPERS_MATH_BLAS_H_
#define SRC_MATH_WRAPPERS_MATH_BLAS_H_

#include "math/primitives/vector_dense.h"
#include "math/primitives/vector_sparse.h"
#include "math/primitives/matrix_dense.h"
#include "math/primitives/matrix_csr.h"
#include "math/primitives/matrix_ijv.h"
#include "basis/utilities/utils.h"

#include <complex>

namespace espreso {
namespace math {
namespace blas {

    // x = y
    template <typename T>
    void copy(const int size, T *x, const int incX, const T *y, const int incY);

    // x *= alpha
    template <typename T>
    void scale(const int size, const T &alpha, T *x, const int incX);

    // x += alpha * y
    template <typename T>
    void add(const int size, T *x, const int incX, const T &alpha, const T *y, const int incY);

    template <typename T>
    T dot(const int size, const T *x, const int incX, const T *y, const int incY);

    template <typename T>
    utils::remove_complex_t<T> norm(const int size, const T *x, const int incX);

    // y = alpha * A * x + beta * y
    template <typename T, typename I>
    void apply(Vector_Dense<T, I> &y, const T &alpha, const Matrix_Dense<T, I> &a, const T &beta, const Vector_Dense<T, I> &x);

    // y = alpha * At * x + beta * y
    template <typename T, typename I>
    void applyT(Vector_Dense<T, I> &y, const T &alpha, const Matrix_Dense<T, I> &a, const T &beta, const Vector_Dense<T, I> &x);

    template <typename T, typename I>
    void apply_hermitian(Vector_Dense<T, I> &y, const T &alpha, const Matrix_Dense<T, I> &a, const T &beta, const Vector_Dense<T, I> &x);

    template <typename T, typename I>
    void AAt(const Matrix_Dense<T, I> &A, Matrix_Dense<T, I> &AAt, bool trans = false);

    // C = alpha * op(A) * op(B) + beta C
    template <typename T, typename I>
    void multiply(T alpha, const Matrix_Dense<T, I> &A, const Matrix_Dense<T, I> &B, T beta, Matrix_Dense<T, I> &C, bool transA = false, bool transB = false);

    template <typename T, typename I>
    void multiply(T alpha, const Matrix_Dense<T, I> &A, const Vector_Dense<T, I> &B, T beta, Vector_Dense<T, I> &C, bool transA = false);

    template <typename T, typename I>
    void multiply(T alpha, const Vector_Dense<T, I> &A, const Vector_Dense<T, I> &B, T beta, T &out);
}
}
}

#endif /* SRC_MATH_WRAPPERS_MATH_BLAS_H_ */
