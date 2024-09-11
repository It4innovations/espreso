
#include "esinfo/eslog.h"
#include "math/wrappers/math.lapack.h"

#include <complex>

#ifndef HAVE_MKL
#ifdef HAVE_LAPACK
#include "lapacke.h"

namespace espreso {
namespace math {
namespace lapack {

template <>
void solve_sym_upper(Matrix_Dense<double> &A, Matrix_Dense<double> &rhs)
{
    int *ipiv = new int[A.nrows];
    LAPACKE_dsytrf(LAPACK_ROW_MAJOR, 'U', A.nrows, A.vals, A.ncols, ipiv);
    LAPACKE_dsytrs(LAPACK_ROW_MAJOR, 'U', A.nrows, rhs.ncols, A.vals, A.ncols, ipiv, rhs.vals, rhs.ncols);
    delete[] ipiv;
}

template <>
void solve_general(Matrix_Dense<double> &A, Matrix_Dense<double> &rhs)
{
    int *ipiv = new int[A.nrows];
    LAPACKE_dgesv(LAPACK_ROW_MAJOR, A.nrows, rhs.ncols, A.vals, A.ncols, ipiv, rhs.vals, rhs.ncols);
    delete[] ipiv;
}

template <>
void get_eig_sym(Matrix_Dense<double> &A, Vector_Dense<double> &values)
{
    switch (A.shape) {
    case Matrix_Shape::FULL: break;
    case Matrix_Shape::UPPER: LAPACKE_dspev (LAPACK_ROW_MAJOR, 'N', 'U', A.nrows, A.vals, values.vals, nullptr, A.ncols); break;
    case Matrix_Shape::LOWER: LAPACKE_dspev (LAPACK_ROW_MAJOR, 'N', 'L', A.nrows, A.vals, values.vals, nullptr, A.ncols); break;
    }
}

template <>
void get_eig_sym(Matrix_Dense<double> &A, Vector_Dense<double> &values, Matrix_Dense<double> &vectors)
{
    // LAPACK_COL_MAJOR and swap 'L' and 'U' to get vectors in rows
    switch (A.shape) {
    case Matrix_Shape::FULL: break;
    case Matrix_Shape::UPPER: LAPACKE_dspev (LAPACK_COL_MAJOR, 'V', 'L', A.nrows, A.vals, values.vals, vectors.vals, A.ncols); break;
    case Matrix_Shape::LOWER: LAPACKE_dspev (LAPACK_COL_MAJOR, 'V', 'U', A.nrows, A.vals, values.vals, vectors.vals, A.ncols); break;
    }
}

template <>
void get_eig_sym(Matrix_Dense<double> &A, Vector_Dense<double> &values, int begin, int end)
{
    int m; Vector_Dense<int, int> fail; fail.resize(A.ncols);
    switch (A.shape) {
    case Matrix_Shape::FULL: break;
    case Matrix_Shape::UPPER: LAPACKE_dspevx(LAPACK_ROW_MAJOR, 'N', 'I', 'U', A.nrows, A.vals, 0, 0, begin, end, 2 * LAPACKE_dlamch('S'), &m, values.vals, nullptr, A.ncols, fail.vals); break;
    case Matrix_Shape::LOWER: LAPACKE_dspevx(LAPACK_ROW_MAJOR, 'N', 'I', 'L', A.nrows, A.vals, 0, 0, begin, end, 2 * LAPACKE_dlamch('S'), &m, values.vals, nullptr, A.ncols, fail.vals); break;
    }
}

template <>
void get_eig_sym(Matrix_Dense<double> &A, Vector_Dense<double> &values, Matrix_Dense<double> &vectors, int begin, int end)
{
    int m; Vector_Dense<int, int> fail; fail.resize(A.ncols);
    switch (A.shape) {
    case Matrix_Shape::FULL: break;
    case Matrix_Shape::UPPER: LAPACKE_dspevx(LAPACK_ROW_MAJOR, 'V', 'I', 'U', A.nrows, A.vals, 0, 0, begin, end, 2 * LAPACKE_dlamch('S'), &m, values.vals, vectors.vals, A.ncols, fail.vals); break;
    case Matrix_Shape::LOWER: LAPACKE_dspevx(LAPACK_ROW_MAJOR, 'V', 'I', 'L', A.nrows, A.vals, 0, 0, begin, end, 2 * LAPACKE_dlamch('S'), &m, values.vals, vectors.vals, A.ncols, fail.vals); break;
    }
}


}
}
}

#endif
#endif
