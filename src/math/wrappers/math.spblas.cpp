
#include "math/wrappers/math.spblas.h"
#include "math/wrappers/math.spsolver.h"
#include "esinfo/eslog.h"

#include <algorithm>
#include <complex>
#include <type_traits>
#include <vector>

namespace espreso {

#ifndef HAVE_MKL
#ifndef HAVE_SUITESPARSE



template <template <typename, typename, typename> class Matrix, typename T, typename I>
SpBLAS<Matrix, T, I>::SpBLAS()
: matrix{}, _spblas{}
{

}

template <template <typename, typename, typename> class Matrix, typename T, typename I>
SpBLAS<Matrix, T, I>::~SpBLAS()
{

}

template <template <typename, typename, typename> class Matrix, typename T, typename I>
SpBLAS<Matrix, T, I>::SpBLAS(const MatrixType &a)
: matrix{}, _spblas{}
{

}

template <template <typename, typename, typename> class Matrix, typename T, typename I>
void SpBLAS<Matrix, T, I>::insert(const MatrixType &a)
{

}

template <template <typename, typename, typename> class Matrix, typename T, typename I>
void SpBLAS<Matrix, T, I>::insertTransposed(const MatrixType &a)
{

}

template <template <typename, typename, typename> class Matrix, typename T, typename I>
void SpBLAS<Matrix, T, I>::extractUpper(MatrixType &a)
{

}

template <template <typename, typename, typename> class Matrix, typename T, typename I>
void SpBLAS<Matrix, T, I>::apply(Vector_Dense<T> &y, const T &alpha, const T &beta, const Vector_Dense<T> &x)
{
	eslog::error("calling of empty SpBLAS wrapper.\n");
}

template <template <typename, typename, typename> class Matrix, typename T, typename I>
void SpBLAS<Matrix, T, I>::transposeTo(SpBLAS<Matrix, T, I> &A)
{

}

template <template <typename, typename, typename> class Matrix, typename T, typename I>
void SpBLAS<Matrix, T, I>::multiply(SpBLAS<Matrix, T, I> &A, SpBLAS<Matrix, T, I> &B)
{
	eslog::error("calling of empty SpBLAS wrapper.\n");
}

#endif
#endif

template struct SpBLAS<Matrix_CSR, float>;
template struct SpBLAS<Matrix_CSR, double>;
template struct SpBLAS<Matrix_CSR, std::complex<double>>;

}
