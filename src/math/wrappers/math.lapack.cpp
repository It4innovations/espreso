
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


}
}
}

#endif
#endif
