
#include "math.blas.h"
#include "esinfo/eslog.h"

#include <complex>

namespace espreso {
namespace math {
namespace blas {

template <>
int dot(const esint size, const int *x, const esint incX, const int *y, const esint incY)
{
    int sum = 0;
    for (int i = 0; i < size; ++i) {
        sum += x[i * incX] * y[i * incY];
    }
    return sum;
}

template <>
void copy(const esint size, int *x, const esint incX, const int *y, const esint incY)
{
    for (int i = 0; i < size; ++i) {
        x[i * incX] = y[i * incY];
    }
}

template <>
void add(const esint size, int *x, const esint incX, const int &alpha, const int *y, const esint incY)
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
void copy(const esint size, double *x, const esint incX, const double *y, const esint incY)
{
    eslog::error("calling of empty BLAS wrapper.\n");
}

template <>
void copy(const esint size, std::complex<double> *x, const esint incX, const std::complex<double> *y, const esint incY)
{
    eslog::error("calling of empty BLAS wrapper.\n");
}

template <>
void scale(const esint size, const float &alpha, float *x, const esint incX)
{
    eslog::error("calling of empty BLAS wrapper.\n");
}

template <>
void scale(const esint size, const double &alpha, double *x, const esint incX)
{
    eslog::error("calling of empty BLAS wrapper.\n");
}

template <>
void scale(const esint size, const std::complex<double> &alpha, std::complex<double> *x, const esint incX)
{
    eslog::error("calling of empty BLAS wrapper.\n");
}

template <>
void add(const esint size, double *x, const esint incX, const double &alpha, const double *y, const esint incY)
{
    eslog::error("calling of empty BLAS wrapper.\n");
}

template <>
void add(const esint size, std::complex<double> *x, const esint incX, const std::complex<double> &alpha, const std::complex<double> *y, const esint incY)
{
    eslog::error("calling of empty BLAS wrapper.\n");
}

template <>
double dot(const esint size, const double *x, const esint incX, const double *y, const esint incY)
{
    eslog::error("calling of empty BLAS wrapper.\n");
    return 0;
}

template <>
std::complex<double> dot(const esint size, const std::complex<double> *x, const esint incX, const std::complex<double> *y, const esint incY)
{
    eslog::error("calling of empty BLAS wrapper.\n");
    return 0;
}

template <>
float norm(const esint size, const float *x, const esint incX)
{
        eslog::error("calling of empty BLAS wrapper.\n");
        return 0;
}

template <>
double norm(const esint size, const double *x, const esint incX)
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

}
}
}

#endif
#endif

