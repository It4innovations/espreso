
#include "esinfo/eslog.h"
#include "math/wrappers/math.lapack.h"

#include <complex>

#ifndef ESPRESO_USE_WRAPPER_LAPACK_MKL
#ifndef ESPRESO_USE_WRAPPER_LAPACK_LAPACK

namespace espreso {
namespace math {
namespace lapack {

template <>
void solve_sym_upper(Matrix_Dense<double> &A, Matrix_Dense<double> &rhs)
{
    eslog::error("calling of empty LAPACK wrapper.\n");
}

template <>
void solve_general(Matrix_Dense<double> &A, Matrix_Dense<double> &rhs)
{
    eslog::error("calling of empty LAPACK wrapper.\n");
}

template <>
void get_eig_sym(Matrix_Dense<double, int> &A, Vector_Dense<double, int> &values)
{
    eslog::error("calling of empty LAPACK wrapper.\n");
}

template <>
void get_eig_sym(Matrix_Dense<double, int> &A, Vector_Dense<double, int> &values, Matrix_Dense<double, int> &vectors)
{
    eslog::error("calling of empty LAPACK wrapper.\n");
}

template <>
void get_eig_sym(Matrix_Dense<double, int> &A, Vector_Dense<double, int> &values, int begin, int end)
{
    eslog::error("calling of empty LAPACK wrapper.\n");
}

template <>
void get_eig_sym(Matrix_Dense<double, int> &A, Vector_Dense<double, int> &values, Matrix_Dense<double, int> &vectors, int begin, int end)
{
    eslog::error("calling of empty LAPACK wrapper.\n");
}

template <>
void submatrix(const Matrix_Dense<double, int> &input, Matrix_Dense<double, int> &output, int start_row, int end_row, int start_col, int end_col)
{
    eslog::error("calling of empty LAPACK wrapper.\n");
}

}
}
}

#endif
#endif
