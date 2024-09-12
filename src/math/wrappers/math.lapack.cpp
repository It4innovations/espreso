
#include "esinfo/eslog.h"
#include "math/wrappers/math.lapack.h"

#include <complex>

#ifndef HAVE_MKL
#ifndef HAVE_LAPACK

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

}
}
}

#endif
#endif
