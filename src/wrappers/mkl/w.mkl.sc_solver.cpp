#include "math/wrappers/math.sc_solver.h"

#ifdef HAVE_MKL

#include <memory>
#include <mkl.h>
#include "wrappers/pardiso/w.pardiso.type.h"



namespace espreso {

template<typename T, typename I>
struct Schur_Complement_Solver_External_Representation
{
    void * pt[64];
    MKL_INT iparm[64];
    MKL_INT mtype;
    MKL_INT msglvl;
    MKL_INT error;
    const Matrix_CSR<T,I> * matrix;
    const Matrix_CSR<T,I> * As[4];
    Matrix_CSR<T,I> concat_matrix;
    Vector_Dense<I,I> perm;
    I other_size;
    I sc_size;
    Vector_Dense<I,I> map_concat;
    Vector_Dense<T,I> x;
    Vector_Dense<T,I> y;
    int stage;
    bool provided_as_four_small_matrices;
};



template<typename T, typename I>
SchurComplementSolver<T,I>::SchurComplementSolver()
{
    ext = std::make_unique<Schur_Complement_Solver_External_Representation<T,I>>();
    ext->msglvl = 0;
    ext->error = 0;
    ext->stage = 0;

    std::fill_n(ext->pt, 64, nullptr);
    std::fill_n(ext->iparm, 64, 0);
    
    ext->iparm[0] = 1; // I did popullate iparm
    ext->iparm[1] = 2; // metis fill reduction
    ext->iparm[34] = 1; // zero-based indexing
    ext->iparm[35] = 1; // schur complement

    ext->provided_as_four_small_matrices = false;
}



template<typename T, typename I>
SchurComplementSolver<T,I>::SchurComplementSolver(SchurComplementSolver && other) = default;



template<typename T, typename I>
SchurComplementSolver<T,I> & SchurComplementSolver<T,I>::operator=(SchurComplementSolver && other) = default;



template<typename T, typename I>
SchurComplementSolver<T,I>::~SchurComplementSolver()
{
    MKL_INT neg_one = -1;
    pardiso(ext->pt, nullptr, nullptr, nullptr, &neg_one, &ext->matrix->nrows, nullptr, nullptr, nullptr, nullptr, nullptr, ext->iparm, &ext->msglvl, nullptr, nullptr, &ext->error);
    if(ext->error != 0) {
        eslog::error("SchurComplementSolver::factorizeSymbolic(): pardiso error %d\n", (int)ext->error);
    }
}



template<typename T, typename I>
void SchurComplementSolver<T,I>::commitMatrix(const Matrix_CSR<T,I> & A, I sc_size)
{
    if(ext->stage != 0) eslog::error("SchurComplementSolver::commitMatrix(): wrong stage\n");
    if(A.nrows != A.ncols) eslog::error("SchurComplementSolver::commitMatrix(): matrix is not square\n");

    ext->mtype = _pardisoType(A);

    ext->matrix = &A;
    ext->sc_size = sc_size;
    ext->other_size = A.nrows - sc_size;

    ext->perm.resize(A.nrows);
    std::fill_n(ext->perm.vals, ext->other_size, I{0});
    std::fill_n(ext->perm.vals + ext->other_size, ext->sc_size, I{1});

    ext->x.resize(A.nrows);
    ext->y.resize(A.nrows);
    std::fill_n(ext->x.vals, ext->x.size, T{0});

    ext->stage = 1;
}



template<typename T, typename I>
void SchurComplementSolver<T,I>::commitMatrix(const Matrix_CSR<T,I> & A11, const Matrix_CSR<T,I> & A12, const Matrix_CSR<T,I> & A21, const Matrix_CSR<T,I> & A22)
{
    if(ext->stage != 0) eslog::error("SchurComplementSolver::commitMatrix(): wrong stage\n");
    if(A11.nrows != A11.ncols || A22.nrows != A22.ncols) eslog::error("SchurComplementSolver::commitMatrix(): A11 or A22 is not square\n");
    if(A12.nrows != A11.nrows || A12.ncols != A22.ncols || A21.nrows != A22.nrows || A21.ncols != A11.ncols) eslog::error("SchurComplementSolver::commitMatrix(): incompatible matrix sizes\n");

    ext->provided_as_four_small_matrices = true;

    auto & As = ext->As;
    As[0] = &A11;
    As[1] = &A12;
    As[2] = &A21;
    As[3] = &A22;

    ext->concat_matrix.resize(As[0]->nrows + As[3]->nrows, As[0]->ncols + As[3]->ncols, As[0]->nnz + As[1]->nnz + As[2]->nnz + As[3]->nnz);
    ext->map_concat.resize(2 * ext->concat_matrix.nnz);
    math::csrMatrixConcat(ext->concat_matrix, *As[0], *As[1], *As[2], *As[3], ext->map_concat, math::CsrMatrixConcatStage::Pattern);
    ext->concat_matrix.type = Matrix_Type::REAL_SYMMETRIC_INDEFINITE;

    commitMatrix(ext->concat_matrix, A22.nrows);
}



template<typename T, typename I>
void SchurComplementSolver<T,I>::factorizeSymbolic()
{
    if(ext->stage != 1) eslog::error("SchurComplementSolver::factorizeSymbolic(): wrong stage\n");

    MKL_INT phase = 11;
    MKL_INT one = 1;
    pardiso(ext->pt, &one, &one, &ext->mtype, &phase, &ext->matrix->nrows, ext->matrix->vals, ext->matrix->rows, ext->matrix->cols, ext->perm.vals, &one, ext->iparm, &ext->msglvl, nullptr, nullptr, &ext->error);
    if(ext->error != 0) {
        eslog::error("SchurComplementSolver::factorizeSymbolic(): pardiso error %d\n", (int)ext->error);
    }

    ext->stage = 2;
}



template<typename T, typename I>
void SchurComplementSolver<T,I>::updateMatrixValues()
{
    if(ext->stage < 2) eslog::error("SchurComplementSolver::updateMatrixValues(): wrong stage\n");

    if(ext->provided_as_four_small_matrices) {
        auto & As = ext->As;
        math::csrMatrixConcat(ext->concat_matrix, *As[0], *As[1], *As[2], *As[3], ext->map_concat, math::CsrMatrixConcatStage::Values);
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

    Matrix_Dense<T,I> sc_tmp;
    sc_tmp.resize(ext->sc_size, ext->sc_size);

    MKL_INT phase = 22;
    MKL_INT one = 1;
    pardiso(ext->pt, &one, &one, &ext->mtype, &phase, &ext->matrix->nrows, ext->matrix->vals, ext->matrix->rows, ext->matrix->cols, ext->perm.vals, &one, ext->iparm, &ext->msglvl, nullptr, sc_tmp.vals, &ext->error);
    if(ext->error != 0) {
        eslog::error("SchurComplementSolver::factorizeNumericAndGetSc(): pardiso error %d\n", (int)ext->error);
    }

    math::copyMatrixDense(sc, sc_tmp, uplo, alpha);

    ext->stage = 4;
}



template<typename T, typename I>
void SchurComplementSolver<T,I>::solveA11(const Vector_Dense<T,I> & rhs, Vector_Dense<T,I> & sol)
{
    if(ext->stage != 4) eslog::error("SchurComplementSolver::solveA11(): wrong stage\n");
    if(rhs.size != ext->other_size || sol.size != ext->other_size) eslog::error("SchurComplementSolver::solveA11(): vectors have wrong size\n");

    MKL_INT phase;
    MKL_INT one = 1;

    std::copy_n(rhs.vals, ext->other_size, ext->x.vals);

    phase = 331;
    pardiso(ext->pt, &one, &one, &ext->mtype, &phase, &ext->matrix->nrows, ext->matrix->vals, ext->matrix->rows, ext->matrix->cols, nullptr, &one, ext->iparm, &ext->msglvl, ext->x.vals, ext->y.vals, &ext->error);
    if(ext->error != 0) {
        eslog::error("SchurComplementSolver::solveA11(): pardiso error %d\n", (int)ext->error);
    }

    std::fill_n(ext->y.vals + ext->other_size, ext->sc_size, T{0.0});

    phase = 333;
    pardiso(ext->pt, &one, &one, &ext->mtype, &phase, &ext->matrix->nrows, ext->matrix->vals, ext->matrix->rows, ext->matrix->cols, nullptr, &one, ext->iparm, &ext->msglvl, ext->y.vals, ext->x.vals, &ext->error);
    if(ext->error != 0) {
        eslog::error("SchurComplementSolver::solveA11(): pardiso error %d\n", (int)ext->error);
    }

    std::copy_n(ext->x.vals, ext->other_size, sol.vals);
}

}

#include "math/wrappers/math.sc_solver.inst.hpp"

#endif
