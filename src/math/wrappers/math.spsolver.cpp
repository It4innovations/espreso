
#include "math/wrappers/math.spsolver.h"
#include "gpu/gpu_management.h"
#include "esinfo/eslog.h"
#include "gpu/gpu_management.h"

#include <complex>

#ifndef ESPRESO_USE_WRAPPER_SPBLAS_MKL
#ifndef ESPRESO_USE_WRAPPER_SPBLAS_SUITESPARSE

namespace espreso {

template<typename T, typename I>
struct Solver_External_Representation
{
};

template <typename T, typename I>
const char* DirectSparseSolver<T, I>::name()
{
    return "NONE";
}

template <typename T, typename I>
bool DirectSparseSolver<T, I>::provideFactors()
{
        return false;
}

template <typename T, typename I>
bool DirectSparseSolver<T, I>::provideSC()
{
        return false;
}

template <typename T, typename I>
Solver_Factors DirectSparseSolver<T, I>::factorsSymmetry()
{
    return Solver_Factors::NONE;
}

template <typename T, typename I>
DirectSparseSolver<T, I>::DirectSparseSolver()
{
    eslog::error("calling of empty sparse solver wrapper.\n");
}

template <typename T, typename I>
DirectSparseSolver<T, I>::DirectSparseSolver(const Matrix_CSR<T, I> &a)
{
    eslog::error("calling of empty sparse solver wrapper.\n");
}

template <typename T, typename I>
DirectSparseSolver<T, I>::DirectSparseSolver(DirectSparseSolver<T, I> &&other) = default;

template <typename T, typename I>
DirectSparseSolver<T, I> & DirectSparseSolver<T, I>::operator=(DirectSparseSolver<T, I> &&other) = default;

template <typename T, typename I>
DirectSparseSolver<T, I>::~DirectSparseSolver()
{

}

template <typename T, typename I>
void DirectSparseSolver<T, I>::commit(const Matrix_CSR<T,I> &a)
{

}

template <typename T, typename I>
void DirectSparseSolver<T, I>::symbolicFactorization()
{

}

template <typename T, typename I>
void DirectSparseSolver<T, I>::numericalFactorization()
{

}

template <typename T, typename I>
void DirectSparseSolver<T, I>::solve(Vector_Dense<T, I> &rhs, Vector_Dense<T, I> &solution)
{

}

template <typename T, typename I>
void DirectSparseSolver<T, I>::solve(Matrix_Dense<T, I> &rhs, Matrix_Dense<T, I> &solution)
{

}

template <typename T, typename I>
void DirectSparseSolver<T, I>::solveForward (Vector_Dense<T, I> &rhs, Vector_Dense<T, I> &solution)
{

}

template <typename T, typename I>
void DirectSparseSolver<T, I>::solveDiagonal(Vector_Dense<T, I> &rhs, Vector_Dense<T, I> &solution)
{

}

template <typename T, typename I>
void DirectSparseSolver<T, I>::solveBackward(Vector_Dense<T, I> &rhs, Vector_Dense<T, I> &solution)
{

}

template <typename T, typename I>
void DirectSparseSolver<T, I>::solveForward (Matrix_Dense<T, I> &rhs, Matrix_Dense<T, I> &solution)
{

}

template <typename T, typename I>
void DirectSparseSolver<T, I>::solveDiagonal(Matrix_Dense<T, I> &rhs, Matrix_Dense<T, I> &solution)
{

}

template <typename T, typename I>
void DirectSparseSolver<T, I>::solveBackward(Matrix_Dense<T, I> &rhs, Matrix_Dense<T, I> &solution)
{

}

template <typename T, typename I>
I DirectSparseSolver<T, I>::getMatrixSize()
{
    return 0;
}

template <typename T, typename I>
I DirectSparseSolver<T, I>::getMatrixNnz()
{
    return 0;
}

template <typename T, typename I>
I DirectSparseSolver<T, I>::getFactorNnz()
{
    return 0;
}

template <typename T, typename I>
template<typename A>
void DirectSparseSolver<T, I>::getFactorL(Matrix_CSR<T, I, A> &/*L*/, bool /*copyPattern*/, bool /*copyValues*/)
{
    eslog::error("calling of empty sparse solver wrapper.\n");
}

template <typename T, typename I>
template<typename A>
void DirectSparseSolver<T, I>::getFactorU(Matrix_CSR<T, I, A> &/*U*/, bool /*copyPattern*/, bool /*copyValues*/)
{
    eslog::error("calling of empty sparse solver wrapper.\n");
}

template <typename T, typename I>
void DirectSparseSolver<T, I>::getPermutation(Permutation<I> &/*perm*/)
{
    eslog::error("calling of empty sparse solver wrapper.\n");
}

template <typename T, typename I>
void DirectSparseSolver<T, I>::getSC(Matrix_Dense<T,I> &/*sc*/, std::vector<I> &/*indices*/, bool /*symmetric_packed*/)
{
    eslog::error("calling of empty sparse solver wrapper.\n");
}

}

#include "math/wrappers/math.spsolver.inst.hpp"

#endif
#endif
