#include "math/wrappers/math.sc_solver.h"

#ifndef ESPRESO_USE_WRAPPER_SPBLAS_MKL
#ifndef ESPRESO_USE_WRAPPER_SPBLAS_SUITESPARSE

namespace espreso {

template<typename T, typename I>
struct Schur_Complement_Solver_External_Representation {};

template<typename T, typename I>
SchurComplementSolver<T,I>::SchurComplementSolver()
{
    eslog::error("Error: empty SchurComplementSolver implementation\n");
}

template<typename T, typename I>
SchurComplementSolver<T,I>::~SchurComplementSolver()
{
    eslog::error("Error: empty SchurComplementSolver implementation\n");
}

template<typename T, typename I>
SchurComplementSolver<T,I>::SchurComplementSolver(SchurComplementSolver && other) = default;

template<typename T, typename I>
SchurComplementSolver<T,I> & SchurComplementSolver<T,I>::operator=(SchurComplementSolver && other) = default;

template <typename T, typename I>
const char * SchurComplementSolver<T,I>::name()
{
    return "EMPTY";
}

template<typename T, typename I>
void SchurComplementSolver<T,I>::commitMatrix(const Matrix_CSR<T,I> & /*A*/, I /*sc_size*/)
{
    eslog::error("Error: empty SchurComplementSolver implementation\n");
}

template<typename T, typename I>
void SchurComplementSolver<T,I>::commitMatrix(const Matrix_CSR<T,I> & /*A11*/, const Matrix_CSR<T,I> & /*A12*/, const Matrix_CSR<T,I> & /*A21*/, const Matrix_CSR<T,I> & /*A22*/)
{
    eslog::error("Error: empty SchurComplementSolver implementation\n");
}

template<typename T, typename I>
void SchurComplementSolver<T,I>::factorizeSymbolic()
{
    eslog::error("Error: empty SchurComplementSolver implementation\n");
}

template<typename T, typename I>
void SchurComplementSolver<T,I>::updateMatrixValues()
{
    eslog::error("Error: empty SchurComplementSolver implementation\n");
}

template<typename T, typename I>
template<typename A>
void SchurComplementSolver<T,I>::factorizeNumericAndGetSc(Matrix_Dense<T,I,A> & /*sc*/, char /*uplo*/, T /*alpha*/)
{
    eslog::error("Error: empty SchurComplementSolver implementation\n");
}

template<typename T, typename I>
void SchurComplementSolver<T,I>::solveA11(const Vector_Dense<T,I> & /*rhs*/, Vector_Dense<T,I> & /*sol*/)
{
    eslog::error("Error: empty SchurComplementSolver implementation\n");
}

}

#include "math/wrappers/math.sc_solver.inst.hpp"

#endif
#endif
