
#include "math/math.h"
#include "math/primitives/matrix_csr.h"
#include "math/wrappers/math.spblas.h"
#include "esinfo/eslog.h"

#include <complex>

#ifndef ESPRESO_USE_WRAPPER_SPBLAS_MKL
#ifndef ESPRESO_USE_WRAPPER_SPBLAS_SUITESPARSE

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
SpBLAS<Matrix, T, I>::SpBLAS(MatrixType &a, bool trans)
{
    eslog::error("calling of empty sparse blas wrapper.\n");
}

template <template <typename, typename, typename> class Matrix, typename T, typename I>
void SpBLAS<Matrix, T, I>::insert(MatrixType &a, bool trans)
{

}

template <template <typename, typename, typename> class Matrix, typename T, typename I>
void SpBLAS<Matrix, T, I>::apply(Vector_Dense<T, I> &y, const T &alpha, const T &beta, const Vector_Dense<T, I> &x)
{
    eslog::error("call empty SpBLAS wrapper\n");
}

template <template <typename, typename, typename> class Matrix, typename T, typename I>
void SpBLAS<Matrix, T, I>::apply(Matrix_Dense<T, I> &y, const T &alpha, const T &beta, const Matrix_Dense<T, I> &x, bool trans)
{
    eslog::error("call empty SpBLAS wrapper\n");
}



namespace math {
namespace spblas {

struct _handle_trsm {};

template<typename T, typename I>
void trsm(MatrixCsxView_new<T,I> & /*A*/, MatrixDenseView_new<T> & /*X*/, MatrixDenseView_new<T> & /*Y*/, handle_trsm & /*handle*/, char /*stage*/)
{
    eslog::error("call empty SpBLAS wrapper\n");
}

struct _handle_mm {};

template<typename T, typename I>
void mm(MatrixCsxView_new<T,I> & /*A*/, MatrixDenseView_new<T> & /*B*/, MatrixDenseView_new<T> & /*C*/, T /*alpha*/, T /*beta*/, handle_mm & /*handle*/, char /*stage*/)
{
    eslog::error("call empty SpBLAS wrapper\n");
}

}
}

}

#include "math/wrappers/math.spblas.inst.hpp"

#endif
#endif
