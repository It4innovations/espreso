
#include "math/math.h"
#include "math/wrappers/math.solver.h"

#ifdef HAVE_MKL
#ifdef ESPRESO_USE_WRAPPER_DNSOLVER_MKL

#include "mkl_lapacke.h"

namespace espreso {

template <typename T, typename I>
DenseSolver<T, I>::DenseSolver()
{

}

template <typename T, typename I>
DenseSolver<T, I>::DenseSolver(const Matrix_Dense<T, I> &a)
{
    math::copy(this->a, a);
}

template <typename T, typename I>
DenseSolver<T, I>::DenseSolver(Matrix_Dense<T, I> &&a)
: a(std::move(a))
{

}

template <typename T, typename I>
DenseSolver<T, I>::~DenseSolver()
{

}

template <typename T, typename I>
void DenseSolver<T, I>::commit(const Matrix_Dense<T,I> &a)
{
    this->a.resize(a.nrows, a.ncols);
    math::copy(this->a, a);
}

template <typename T, typename I>
void DenseSolver<T, I>::commit(Matrix_Dense<T,I> &&a)
{
    this->a = std::move(a);
}

template <typename T, typename I>
void DenseSolver<T, I>::factorization()
{
    ipiv.resize(a.nrows);
    LAPACKE_dsytrf(LAPACK_ROW_MAJOR, 'U', a.nrows, a.vals, a.ncols, ipiv.data());
}

template <typename T, typename I>
void DenseSolver<T, I>::solve(Vector_Dense<T, I> &rhs)
{
    LAPACKE_dsytrs(LAPACK_ROW_MAJOR, 'U', a.nrows, 1, a.vals, a.ncols, ipiv.data(), rhs.vals, 1);
}

template <typename T, typename I>
void DenseSolver<T, I>::solve(Matrix_Dense<T, I> &rhs)
{
    LAPACKE_dsytrs(LAPACK_ROW_MAJOR, 'U', a.nrows, rhs.ncols, a.vals, a.ncols, ipiv.data(), rhs.vals, rhs.ncols);
}

template struct DenseSolver<double, int>;

}

#endif
#endif
