
#include "math.solver.h"

#ifndef HAVE_MKL
#ifndef HAVE_LAPACK

namespace espreso {

template <typename T, typename I>
DenseSolver<T, I>::DenseSolver()
{
    eslog::error("calling of empty sparse solver wrapper.\n");
}

template <typename T, typename I>
DenseSolver<T, I>::DenseSolver(const Matrix_Dense<T, I> &a)
{
    eslog::error("calling of empty sparse solver wrapper.\n");
}

template <typename T, typename I>
DenseSolver<T, I>::DenseSolver(Matrix_Dense<T, I> &&a)
{
    eslog::error("calling of empty sparse solver wrapper.\n");
}

template <typename T, typename I>
DenseSolver<T, I>::~DenseSolver()
{

}

template <typename T, typename I>
void DenseSolver<T, I>::commit(const Matrix_Dense<T,I> &a)
{

}

template <typename T, typename I>
void DenseSolver<T, I>::commit(Matrix_Dense<T,I> &&a)
{

}

template <typename T, typename I>
void DenseSolver<T, I>::factorization()
{

}

template <typename T, typename I>
void DenseSolver<T, I>::solve(Vector_Dense<T, I> &rhs, Vector_Dense<T, I> &solution)
{

}

template <typename T, typename I>
void DenseSolver<T, I>::solve(Matrix_Dense<T, I> &rhs, Matrix_Dense<T, I> &solution)
{

}

template struct DenseSolver<double, int>;

}

#endif
#endif
