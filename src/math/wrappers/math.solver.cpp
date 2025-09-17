
#include "math.solver.h"

#ifndef ESPRESO_USE_WRAPPER_LAPACK_MKL
#ifndef ESPRESO_USE_WRAPPER_LAPACK_LAPACK

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
void DenseSolver<T, I>::solve(Vector_Dense<T, I> &rhs)
{

}

template <typename T, typename I>
void DenseSolver<T, I>::solve(Matrix_Dense<T, I> &rhs)
{

}

}

#include "math/wrappers/math.solver.inst.hpp"

#endif
#endif
