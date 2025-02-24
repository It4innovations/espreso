
#include "math.blas.h"
#include "esinfo/eslog.h"

#include <complex>

namespace espreso {
namespace math {
namespace blas {

template <>
int dot(const int size, const int *x, const int incX, const int *y, const int incY)
{
    int sum = 0;
    for (int i = 0; i < size; ++i) {
        sum += x[i * incX] * y[i * incY];
    }
    return sum;
}

template <>
void copy(const int size, int *x, const int incX, const int *y, const int incY)
{
    for (int i = 0; i < size; ++i) {
        x[i * incX] = y[i * incY];
    }
}

template <>
void add(const int size, int *x, const int incX, const int &alpha, const int *y, const int incY)
{
    for (int i = 0; i < size; ++i) {
        x[i * incX] += alpha * y[i * incY];
    }
}

}
}
}

#ifndef ESPRESO_USE_WRAPPER_DNBLAS_MKL
#ifndef ESPRESO_USE_WRAPPER_DNBLAS_BLAS

namespace espreso {
namespace math {
namespace blas {

template <typename T>
void copy(const int size, T *x, const int incX, const T *y, const int incY)
{
    eslog::error("calling of empty BLAS wrapper.\n");
}

template <typename T>
void scale(const int size, const T &alpha, T *x, const int incX)
{
    eslog::error("calling of empty BLAS wrapper.\n");
}

template <typename T>
void add(const int size, T *x, const int incX, const T &alpha, const T *y, const int incY)
{
    eslog::error("calling of empty BLAS wrapper.\n");
}

template <typename T>
T dot(const int size, const T *x, const int incX, const T *y, const int incY)
{
    eslog::error("calling of empty BLAS wrapper.\n");
    return 0;
}

template <typename T>
utils::remove_complex_t<T> norm(const int size, const T *x, const int incX)
{
    eslog::error("calling of empty BLAS wrapper.\n");
    return 0;
}

template <typename T, typename I>
void apply(Vector_Dense<T, I> &y, const T &alpha, const Matrix_Dense<T, I> &a, const T &beta, const Vector_Dense<T, I> &x)
{
    eslog::error("calling of empty BLAS wrapper.\n");
}

template <typename T, typename I>
void applyT(Vector_Dense<T, I> &y, const T &alpha, const Matrix_Dense<T, I> &a, const T &beta, const Vector_Dense<T, I> &x)
{
    eslog::error("calling of empty BLAS wrapper.\n");
}

template <typename T, typename I>
void apply_hermitian(Vector_Dense<T, I> &y, const T &alpha, const Matrix_Dense<T, I> &a, char uplo, const T &beta, const Vector_Dense<T, I> &x)
{
    eslog::error("calling of empty BLAS wrapper.\n");
}

template <typename T, typename I>
void AAt(const Matrix_Dense<T, I> &A, Matrix_Dense<T, I> &AAt, bool trans)
{
    eslog::error("calling of empty BLAS wrapper.\n");
}

template <typename T, typename I>
void multiply(T alpha, const Matrix_Dense<T, I> &A, const Matrix_Dense<T, I> &B, T beta, Matrix_Dense<T, I> &C, bool transA, bool transB)
{
    eslog::error("calling of empty BLAS wrapper.\n");
}

template <typename T, typename I>
void multiply(T alpha, const Matrix_Dense<T, I> &A, const Vector_Dense<T, I> &B, T beta, Vector_Dense<T, I> &C, bool transA)
{
    eslog::error("calling of empty BLAS wrapper.\n");
}

template <typename T, typename I>
void multiply(T alpha, const Vector_Dense<T, I> &A, const Vector_Dense<T, I> &B, T beta, T &out)
{
    eslog::error("calling of empty BLAS wrapper.\n");
}

template<typename T>
void trsm(MatrixDenseView_new<T> & A, MatrixDenseView_new<T> & X, T alpha)
{
    eslog::error("calling of empty BLAS wrapper.\n");
}

template<typename T>
void gemm(MatrixDenseView_new<T> & A, MatrixDenseView_new<T> & B, MatrixDenseView_new<T> & C, T alpha, T beta)
{
    eslog::error("calling of empty BLAS wrapper.\n");
}

template<typename T>
void herk(MatrixDenseView_new<T> & A, MatrixDenseView_new<T> & C, herk_mode mode, T alpha, T beta)
{
    eslog::error("calling of empty BLAS wrapper.\n");
}

}
}
}

#include "math.blas.inst.hpp"

#endif
#endif
