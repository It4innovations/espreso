#include "math/wrappers/math.sc_solver.h"

#ifdef HAVE_SUITESPARSE

#include "w.suitesparse.cholmod.h"



namespace espreso {

template<typename T, typename I>
static void extract_submatrices(const Matrix_CSR<T,I> & A, Matrix_CSR<T,I> & A11, Matrix_CSR<T,I> & A12, Matrix_CSR<T,I> & A22, I sc_size)
{
    I other_size = A.nrows - sc_size;

    SpBLAS<Matrix_CSR,T,I>::submatrix(A, A11, 0,          other_size, 0,          other_size);
    SpBLAS<Matrix_CSR,T,I>::submatrix(A, A12, 0,          other_size, other_size, A.nrows);
    SpBLAS<Matrix_CSR,T,I>::submatrix(A, A22, other_size, A.nrows,    other_size, A.nrows);
}



template<typename T, typename I>
struct Schur_Complement_Solver_External_Representation
{
    cholmod_common cm_common;
    int stage;
    const Matrix_CSR<T,I> * A11;
    const Matrix_CSR<T,I> * A12;
    const Matrix_CSR<T,I> * A22;
    Matrix_CSR<T,I> * A11_my = nullptr;
    Matrix_CSR<T,I> * A12_my = nullptr;
    Matrix_CSR<T,I> * A22_my = nullptr;
    const Matrix_CSR<T,I> * A;
    Matrix_Dense<T,I> X;
    Matrix_Dense<T,I> sc;
    cholmod_sparse cm_A11_view;
    cholmod_sparse cm_A12_view;
    cholmod_factor * cm_factor;
    cholmod_dense cm_X_view;
    cholmod_dense cm_sc_view;
    I other_size;
    I sc_size;
    bool provided_as_one_large_matrix;
};



template<typename T, typename I>
SchurComplementSolver<T,I>::SchurComplementSolver()
{
    ext = std::make_unique<Schur_Complement_Solver_External_Representation<T,I>>();

    ext->provided_as_one_large_matrix = false;

    _start<I>(ext->cm_common);
    ext->cm_common.final_ll = 1;
    ext->cm_common.nthreads_max = 1;
    ext->cm_common.nmethods = 1;
    ext->cm_common.method[0].ordering = CHOLMOD_METIS;
    ext->cm_common.itype = _getCholmodItype<I>();
    ext->cm_common.supernodal = CHOLMOD_SUPERNODAL;
}



template<typename T, typename I>
SchurComplementSolver<T,I>::~SchurComplementSolver()
{
    _free<I>(ext->cm_factor, ext->cm_common);
    _finish<I>(ext->cm_common);

    delete ext->A11_my;
    delete ext->A12_my;
    delete ext->A22_my;
}



template<typename T, typename I>
SchurComplementSolver<T,I>::SchurComplementSolver(SchurComplementSolver && other) = default;



template<typename T, typename I>
SchurComplementSolver<T,I> & SchurComplementSolver<T,I>::operator=(SchurComplementSolver && other) = default;



template<typename T, typename I>
void SchurComplementSolver<T,I>::commitMatrix(const Matrix_CSR<T,I> & A, I sc_size)
{
    if(ext->stage != 0) eslog::error("SchurComplementSolver::commitMatrix(): wrong stage\n");
    if(A.nrows != A.ncols) eslog::error("SchurComplementSolver::commitMatrix(): matrix is not square\n");
    if(getSymmetry(A.type) != Matrix_Symmetry::HERMITIAN || A.shape != Matrix_Shape::UPPER) eslog::error("SchurComplementSolver::commitMatrix(): only works with hermitian matrices stored in upper triangle now. TODO.\n");

    ext->provided_as_one_large_matrix = true;
    ext->sc_size = sc_size;
    ext->other_size = A.nrows - sc_size;

    ext->A = &A;
    ext->A11_my = new Matrix_CSR<T,I>();
    ext->A12_my = new Matrix_CSR<T,I>();
    ext->A22_my = new Matrix_CSR<T,I>();
    ext->A11_my->shape = Matrix_Shape::UPPER;
    ext->A22_my->shape = Matrix_Shape::UPPER;
    ext->A11_my->type = A.type;
    ext->A22_my->type = A.type;
    ext->A11 = ext->A11_my;
    ext->A12 = ext->A12_my;
    ext->A22 = ext->A22_my;

    extract_submatrices(A, *ext->A11_my, *ext->A12_my, *ext->A22_my, ext->sc_size);
    
    ext->stage = 1;
}



template<typename T, typename I>
void SchurComplementSolver<T,I>::commitMatrix(const Matrix_CSR<T,I> & A11, const Matrix_CSR<T,I> & A12, const Matrix_CSR<T,I> & A21, const Matrix_CSR<T,I> & A22)
{
    if(ext->stage != 0) eslog::error("SchurComplementSolver::commitMatrix(): wrong stage\n");
    if(A11.nrows != A11.ncols || A22.nrows != A22.ncols) eslog::error("SchurComplementSolver::commitMatrix(): A11 or A22 is not square\n");
    if(A12.nrows != A11.nrows || A12.ncols != A22.ncols || A21.nrows != A22.nrows || A21.ncols != A11.ncols) eslog::error("SchurComplementSolver::commitMatrix(): incompatible matrix sizes\n");

    if(getSymmetry(A11.type) != Matrix_Symmetry::HERMITIAN || getSymmetry(A22.type) != Matrix_Symmetry::HERMITIAN) eslog::error("SchurComplementSolver::commitMatrix(): only works with hermitian matrices stored in upper triangle now. TODO.\n");
    if(A11.shape != Matrix_Shape::UPPER || A22.shape != Matrix_Shape::UPPER || A21.nnz != 0) eslog::error("SchurComplementSolver::commitMatrix(): only works with hermitian matrices stored in upper triangle now. TODO.\n");
    if(utils::is_real<T>() && A11.type != Matrix_Type::REAL_SYMMETRIC_POSITIVE_DEFINITE) eslog::error("SchurComplementSolver::commitMatrix(): A11 has to be symmetric positive definite\n");
    if(utils::is_complex<T>() && A11.type != Matrix_Type::COMPLEX_HERMITIAN_POSITIVE_DEFINITE) eslog::error("SchurComplementSolver::commitMatrix(): A11 has to be hermitian positive definite\n");

    ext->sc_size = A22.nrows;
    ext->other_size = A11.nrows;

    ext->A11 = &A11;
    ext->A12 = &A12;
    // ext->A21 = &A21;
    ext->A22 = &A22;

    ext->stage = 1;
}



template<typename T, typename I>
void SchurComplementSolver<T,I>::factorizeSymbolic()
{
    if(ext->stage != 1) eslog::error("SchurComplementSolver::factorizeSymbolic(): wrong stage\n");

    // finish what I did not do in commit
    {
        ext->X.resize(ext->sc_size, ext->other_size);
        math::csrToDense(ext->X, *ext->A12, true, math::CsrToDenseStage::ZeroFill);

        popullateCholmodSparse(ext->cm_A11_view, *ext->A11);
        popullateCholmodSparse(ext->cm_A12_view, *ext->A12);
        popullateCholmodDense(ext->cm_X_view, ext->X);
        ext->sc.resize(ext->sc_size, ext->sc_size);
        popullateCholmodDense(ext->cm_sc_view, ext->sc);
    }

    ext->cm_factor = _analyze<I>(&ext->cm_A11_view, ext->cm_common);

    ext->stage = 2;
}



template<typename T, typename I>
void SchurComplementSolver<T,I>::updateMatrixValues()
{
    if(ext->stage < 2) eslog::error("SchurComplementSolver::updateMatrixValues(): wrong stage\n");

    if(ext->provided_as_one_large_matrix) {
        extract_submatrices(*ext->A, *ext->A11_my, *ext->A12_my, *ext->A22_my, ext->sc_size);
    }

    ext->stage = 3;
}



template<typename T, typename I>
template<typename A>
void SchurComplementSolver<T,I>::factorizeNumericAndGetSc(Matrix_Dense<T,I,A> & sc, char uplo, T alpha)
{
    if(ext->stage != 3) eslog::error("SchurComplementSolver::factorizeNumericAndGetSc(): wrong stage\n");
    if(sc.nrows != ext->sc_size || sc.ncols != ext->sc_size) eslog::error("SchurComplementSolver::factorizeNumericAndGetSc(): output SC matrix has wrong size\n");
    if(sc.shape != Matrix_Shape::FULL) eslog::error("SchurComplementSolver::factorizeNumericAndGetSc(): output SC matrix needs to be allocated as full matrix\n");

    _factorize<I>(ext->cm_factor, &ext->cm_A11_view, ext->cm_common);

    math::csrToDense(ext->X, *ext->A12, true, math::CsrToDenseStage::ValuesCopy);

    cholmod_dense * cm_Y = _solve<I>(CHOLMOD_A, ext->cm_factor, &ext->cm_X_view, ext->cm_common);

    T one[2] = {1.0, 0.0};
    T zero[2] = {0.0, 0.0};
    cholmod_sdmult(&ext->cm_A12_view, 0, one, zero, cm_Y, &ext->cm_sc_view, &ext->cm_common);
    
    math::copyMatrixDense(sc, ext->sc, uplo, alpha);

    _free<I>(cm_Y, ext->cm_common);
    
    ext->stage = 4;
}



template<typename T, typename I>
void SchurComplementSolver<T,I>::solveA11(const Vector_Dense<T,I> & rhs, Vector_Dense<T,I> & sol)
{
    if(ext->stage != 4) eslog::error("SchurComplementSolver::factorizeNumericAndGetSc(): wrong stage\n");
    if(rhs.size != ext->other_size || sol.size != ext->other_size) eslog::error("SchurComplementSolver::solveA11(): vectors have wrong size\n");

    cholmod_dense cm_rhs;
    popullateCholmodDense(cm_rhs, rhs);

    cholmod_dense * cm_sol = _solve<I>(CHOLMOD_A, ext->cm_factor, &cm_rhs, ext->cm_common);

    std::copy_n((T*)cm_sol->x, sol.size, sol.vals);
    
    _free<I>(cm_sol, ext->cm_common);
}



}

#include "math/wrappers/math.sc_solver.inst.hpp"

#endif
