
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

#ifndef HAVE_MKL
#ifndef HAVE_BLAS

namespace espreso {
namespace math {
namespace blas {

template <>
void copy(const int size, double *x, const int incX, const double *y, const int incY)
{
    eslog::error("calling of empty BLAS wrapper.\n");
}

template <>
void copy(const int size, std::complex<double> *x, const int incX, const std::complex<double> *y, const int incY)
{
    eslog::error("calling of empty BLAS wrapper.\n");
}

template <>
void scale(const int size, const float &alpha, float *x, const int incX)
{
    eslog::error("calling of empty BLAS wrapper.\n");
}

template <>
void scale(const int size, const double &alpha, double *x, const int incX)
{
    eslog::error("calling of empty BLAS wrapper.\n");
}

template <>
void scale(const int size, const std::complex<double> &alpha, std::complex<double> *x, const int incX)
{
    eslog::error("calling of empty BLAS wrapper.\n");
}

template <>
void add(const int size, double *x, const int incX, const double &alpha, const double *y, const int incY)
{
    eslog::error("calling of empty BLAS wrapper.\n");
}

template <>
void add(const int size, std::complex<double> *x, const int incX, const std::complex<double> &alpha, const std::complex<double> *y, const int incY)
{
    eslog::error("calling of empty BLAS wrapper.\n");
}

template <>
double dot(const int size, const double *x, const int incX, const double *y, const int incY)
{
    eslog::error("calling of empty BLAS wrapper.\n");
    return 0;
}

template <>
std::complex<double> dot(const int size, const std::complex<double> *x, const int incX, const std::complex<double> *y, const int incY)
{
    eslog::error("calling of empty BLAS wrapper.\n");
    return 0;
}

template <>
float norm(const int size, const float *x, const int incX)
{
        eslog::error("calling of empty BLAS wrapper.\n");
        return 0;
}

template <>
double norm(const int size, const double *x, const int incX)
{
    eslog::error("calling of empty BLAS wrapper.\n");
    return 0;
}

template <>
float norm(const int size, const std::complex<float> *x, const int incX)
{
    eslog::error("calling of empty BLAS wrapper.\n");
    return 0;
}

template <>
double norm(const int size, const std::complex<double> *x, const int incX)
{
    eslog::error("calling of empty BLAS wrapper.\n");
    return 0;
}

template <>
void apply(Vector_Dense<double, int> &y, const double &alpha, const Matrix_Dense<double, int> &a, const double &beta, const Vector_Dense<double, int> &x)
{
    eslog::error("calling of empty BLAS wrapper.\n");
}

template <>
void apply(Vector_Dense<std::complex<double>, int> &y, const std::complex<double> &alpha, const Matrix_Dense<std::complex<double>, int> &a, const std::complex<double> &beta, const Vector_Dense<std::complex<double>, int> &x)
{
    eslog::error("calling of empty BLAS wrapper.\n");
}

template <>
void applyT(Vector_Dense<double, int> &y, const double &alpha, const Matrix_Dense<double, int> &a, const double &beta, const Vector_Dense<double, int> &x)
{
    eslog::error("calling of empty BLAS wrapper.\n");
}

template <>
void applyT(Vector_Dense<std::complex<double> > &y, const std::complex<double> &alpha, const Matrix_Dense<std::complex<double> > &a, const std::complex<double> &beta, const Vector_Dense<std::complex<double> > &x)
{
    eslog::error("calling of empty BLAS wrapper.\n");
}

template <>
void AAt(const Matrix_Dense<double> &A, Matrix_Dense<double> &AAt, bool trans)
{
    eslog::error("calling of empty BLAS wrapper.\n");
}

template <>
void multiply(double alpha, const Matrix_Dense<double> &A, const Matrix_Dense<double> &B, double beta, Matrix_Dense<double> &C, bool transA, bool transB)
{
    eslog::error("calling of empty BLAS wrapper.\n");
}

template <>
void multiply(double alpha, const Matrix_Dense<double> &A, const Vector_Dense<double> &B, double beta, Vector_Dense<double> &C, bool transA)
{
    eslog::error("calling of empty BLAS wrapper.\n");
}

}
}
}

#endif
#endif

