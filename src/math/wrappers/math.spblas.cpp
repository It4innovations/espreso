
#include "math/wrappers/math.spblas.h"
#include "esinfo/eslog.h"

#include <complex>

#ifndef HAVE_MKL
#ifndef HAVE_SUITESPARSE

namespace espreso {

struct Matrix_SpBLAS_External_Representation {

};

template <template <typename, typename, typename> class Matrix, typename T, typename I>
SpBLAS<Matrix, T, I>::SpBLAS()
{
    eslog::error("calling of empty sparse blas wrapper.\n");
}

template <template <typename, typename, typename> class Matrix, typename T, typename I>
SpBLAS<Matrix, T, I>::~SpBLAS()
{

}

template <template <typename, typename, typename> class Matrix, typename T, typename I>
SpBLAS<Matrix, T, I>::SpBLAS(MatrixType &a)

{

}

template <template <typename, typename, typename> class Matrix, typename T, typename I>
void SpBLAS<Matrix, T, I>::insert(MatrixType &a)
{

}

template <>
void SpBLAS<Matrix_CSR, float, int>::apply(Vector_Dense<float> &y, const float &alpha, const float &beta, const Vector_Dense<float> &x)
{

}

template <>
void SpBLAS<Matrix_CSR, double, int>::apply(Vector_Dense<double> &y, const double &alpha, const double &beta, const Vector_Dense<double> &x)
{

}

template <>
void SpBLAS<Matrix_CSR, std::complex<double>, int>::apply(Vector_Dense<std::complex<double>, int> &y, const std::complex<double> &alpha, const std::complex<double> &beta, const Vector_Dense<std::complex<double> > &x)
{

}

template struct SpBLAS<Matrix_CSR, float, int>;
template struct SpBLAS<Matrix_CSR, double, int>;
template struct SpBLAS<Matrix_CSR, std::complex<double>, int>;

}

#endif
#endif
