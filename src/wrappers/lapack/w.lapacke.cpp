
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


}
}
}

#endif
#endif
