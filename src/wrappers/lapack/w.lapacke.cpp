
#include "math/math.h"
#include "esinfo/eslog.h"

#include <complex>

#ifdef HAVE_LAPACK
#ifdef ESPRESO_USE_WRAPPER_LAPACK_LAPACK

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
    int ret = 0;
    switch (A.shape) {
    case Matrix_Shape::FULL:  ret = LAPACKE_dsyev(LAPACK_ROW_MAJOR, 'V', 'U', A.nrows, A.vals, A.ncols, values.vals); break;
    case Matrix_Shape::UPPER: ret = LAPACKE_dspev(LAPACK_ROW_MAJOR, 'N', 'U', A.nrows, A.vals, values.vals, nullptr, A.ncols); break;
    case Matrix_Shape::LOWER: ret = LAPACKE_dspev(LAPACK_ROW_MAJOR, 'N', 'L', A.nrows, A.vals, values.vals, nullptr, A.ncols); break;
    }
}

template <>
void get_eig_sym(Matrix_Dense<double> &A, Vector_Dense<double> &values, Matrix_Dense<double> &vectors)
{
    int ret = 0;
    // LAPACK_COL_MAJOR and swap 'L' and 'U' to get vectors in rows
    switch (A.shape) {
    case Matrix_Shape::FULL: math::copy(vectors, A); ret = LAPACKE_dsyev(LAPACK_ROW_MAJOR, 'V', 'U', A.nrows, vectors.vals, A.ncols, values.vals); break;
    case Matrix_Shape::UPPER: ret = LAPACKE_dspev(LAPACK_ROW_MAJOR, 'V', 'U', A.nrows, A.vals, values.vals, vectors.vals, A.ncols); break;
    case Matrix_Shape::LOWER: ret = LAPACKE_dspev(LAPACK_ROW_MAJOR, 'V', 'L', A.nrows, A.vals, values.vals, vectors.vals, A.ncols); break;
    }
}

template <>
void get_eig_sym(Matrix_Dense<double> &A, Vector_Dense<double> &values, int begin, int end)
{
    int m; Vector_Dense<int, int> fail; fail.resize(2 * A.ncols);
    int ret = 0;
    switch (A.shape) {
    case Matrix_Shape::FULL:  ret = LAPACKE_dsyevr(LAPACK_ROW_MAJOR, 'N', 'I', 'U', A.nrows, A.vals, A.ncols, 0, 0, begin, end, LAPACKE_dlamch('S'), &m, values.vals, nullptr, A.ncols, fail.vals); break;
    case Matrix_Shape::UPPER: ret = LAPACKE_dspevx(LAPACK_ROW_MAJOR, 'N', 'I', 'U', A.nrows, A.vals, 0, 0, begin, end, 2 * LAPACKE_dlamch('S'), &m, values.vals, nullptr, A.ncols, fail.vals); break;
    case Matrix_Shape::LOWER: ret = LAPACKE_dspevx(LAPACK_ROW_MAJOR, 'N', 'I', 'L', A.nrows, A.vals, 0, 0, begin, end, 2 * LAPACKE_dlamch('S'), &m, values.vals, nullptr, A.ncols, fail.vals); break;
    }
}

template <>
void get_eig_sym(Matrix_Dense<double> &A, Vector_Dense<double> &values, Matrix_Dense<double> &vectors, int begin, int end)
{
    int m; Vector_Dense<int, int> fail; fail.resize(2 * A.ncols);
    int ret = 0;
    switch (A.shape) {
    case Matrix_Shape::FULL:  ret = LAPACKE_dsyevr(LAPACK_ROW_MAJOR, 'V', 'I', 'U', A.nrows, A.vals, A.ncols, 0, 0, begin, end, LAPACKE_dlamch('S'), &m, values.vals, vectors.vals, vectors.ncols, fail.vals); break;
    case Matrix_Shape::UPPER: ret = LAPACKE_dspevx(LAPACK_ROW_MAJOR, 'V', 'I', 'U', A.nrows, A.vals, 0, 0, begin, end, 2 * LAPACKE_dlamch('S'), &m, values.vals, vectors.vals, A.ncols, fail.vals); break;
    case Matrix_Shape::LOWER: ret = LAPACKE_dspevx(LAPACK_ROW_MAJOR, 'V', 'I', 'L', A.nrows, A.vals, 0, 0, begin, end, 2 * LAPACKE_dlamch('S'), &m, values.vals, vectors.vals, A.ncols, fail.vals); break;
    }
}

template <>
void submatrix(const Matrix_Dense<double, int> &input, Matrix_Dense<double, int> &output, int start_row, int end_row, int start_col, int end_col)
{
    if (input.shape != Matrix_Shape::FULL) eslog::error("Cannot copy triangular matrix.\n");

    output.resize(end_row - start_row, end_col - start_col);
    LAPACKE_dlacpy(LAPACK_ROW_MAJOR, ' ', end_row - start_row, end_col - start_col, input.vals + start_row * input.ncols + start_col, input.ncols, output.vals, output.ncols);
}


}
}
}

#endif
#endif
