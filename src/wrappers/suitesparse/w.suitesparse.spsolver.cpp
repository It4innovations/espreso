
#include "math/math.h"
#include "esinfo/eslog.h"
#include "gpu/gpu_management.h"
#include "basis/utilities/utils.h"

#include <complex>

#ifdef HAVE_SUITESPARSE

#include "math/wrappers/math.spsolver.h"
#include "w.suitesparse.cholmod.h"

namespace espreso {

template<typename T, typename I>
struct Solver_External_Representation {
    cholmod_common cm_common;
    cholmod_factor * cm_factor_super = nullptr;
    cholmod_factor * cm_factor_simpl = nullptr;
    cholmod_sparse * cm_matrix_view = nullptr;
    const Matrix_CSR<T, I> * matrix = nullptr;
    Vector_Dense<I> map_simpl_super;
    I factor_nnz;
    int stage = 0; // 0 = completely uninitialized, 1 = initialized without matrix, 2 = have matrix, 3 = symbolic factorization done, 4 = numeric factorization done
    bool getfactor_preprocess_done = false;
    bool getfactor_preprocess_lazy;
    char zerodrop;
};

template<typename T, typename I>
void getfactor_preprocess(std::unique_ptr<Solver_External_Representation<T,I>> & ext)
{
    if(ext->cm_factor_super->xsize > utils::get_max_val_no_precision_loss_in_fp<T>()) eslog::error("symbolicFactorization: factor nnz too large for my super->simpl map\n");

    ext->cm_factor_simpl = _copyFactor<I>(ext->cm_factor_super, ext->cm_common);

    _changeFactor<I>(_getCholmodXtype<T>(), true, true, true, true, ext->cm_factor_simpl, ext->cm_common);

    for(size_t i = 0; i < ext->cm_factor_simpl->xsize; i++) {
        reinterpret_cast<T*>(ext->cm_factor_simpl->x)[i] = static_cast<T>(i);
    }

    _changeFactor<I>(_getCholmodXtype<T>(), true, false, true, true, ext->cm_factor_simpl, ext->cm_common);

    if(ext->zerodrop == 'D') {
        _resymbol<I>(ext->cm_matrix_view, nullptr, 0, 1, ext->cm_factor_simpl, ext->cm_common);
    }

    if(ext->cm_factor_simpl->nzmax != (size_t)ext->factor_nnz) {
        eslog::error("getfactor_preprocess: some weird error in cholmod\n");
    }

    ext->map_simpl_super.resize(ext->factor_nnz);
    for(I i = 0; i < ext->map_simpl_super.size; i++) {
        ext->map_simpl_super.vals[i] = static_cast<I>(std::real(reinterpret_cast<T*>(ext->cm_factor_simpl->x)[i]));
    }

    ext->getfactor_preprocess_done = true;
}



template <typename T, typename I>
const char* DirectSparseSolver<T, I>::name()
{
    return "SUITESPARSE";
}

template <typename T, typename I>
bool DirectSparseSolver<T, I>::provideFactors()
{
    return true;
}

template <typename T, typename I>
bool DirectSparseSolver<T, I>::provideSC()
{
    // manually computed
    return true;
}

template <typename T, typename I>
Solver_Factors DirectSparseSolver<T, I>::factorsSymmetry()
{
    return Solver_Factors::HERMITIAN_UPPER;
}

template <typename T, typename I>
DirectSparseSolver<T, I>::DirectSparseSolver()
{
    ext = std::make_unique<Solver_External_Representation<T,I>>();

    ext->zerodrop = 'D'; // Drop, Keep
    ext->getfactor_preprocess_lazy = true;

    _start<I>(ext->cm_common);
    ext->cm_common.final_ll = 1;
    ext->cm_common.nthreads_max = 1;
    ext->cm_common.nmethods = 1;
    ext->cm_common.method[0].ordering = CHOLMOD_METIS;
    ext->cm_common.itype = _getCholmodItype<I>();
    ext->cm_common.supernodal = CHOLMOD_SUPERNODAL;

    ext->stage = 1;
}

template <typename T, typename I>
DirectSparseSolver<T, I>::DirectSparseSolver(const Matrix_CSR<T, I> &a) : DirectSparseSolver()
{
    commit(a);
}

template <typename T, typename I>
DirectSparseSolver<T, I>::DirectSparseSolver(DirectSparseSolver<T, I> &&other) = default;

template <typename T, typename I>
DirectSparseSolver<T, I> & DirectSparseSolver<T, I>::operator=(DirectSparseSolver<T, I> &&other) = default;

template <typename T, typename I>
DirectSparseSolver<T, I>::~DirectSparseSolver()
{
    delete ext->cm_matrix_view;

    if(ext->cm_factor_super != nullptr) _free<I>(ext->cm_factor_super, ext->cm_common);

    if(ext->cm_factor_simpl != nullptr) _free<I>(ext->cm_factor_simpl, ext->cm_common);

    _finish<I>(ext->cm_common);
}

template <typename T, typename I>
void DirectSparseSolver<T, I>::commit(const Matrix_CSR<T,I> &a)
{
    if(a.nrows != a.ncols) eslog::error("commit: matrix has to be square\n");
    if(a.type != Matrix_Type::REAL_SYMMETRIC_POSITIVE_DEFINITE && a.type != Matrix_Type::COMPLEX_HERMITIAN_POSITIVE_DEFINITE) eslog::error("commit: matrix has to be SPD or HPD\n");
    if(a.shape != Matrix_Shape::UPPER) eslog::error("commit: CSR matrix has to be upper triangular\n");

    if(ext->stage < 1) eslog::error("commit: invalid order of operations in spsolver\n");

    ext->matrix = &a;

    if(ext->cm_matrix_view == nullptr) {
        ext->cm_matrix_view = new cholmod_sparse();
        ext->cm_matrix_view->nrow = a.ncols;
        ext->cm_matrix_view->ncol = a.nrows;
        ext->cm_matrix_view->nzmax = a.nnz;
        ext->cm_matrix_view->nz = nullptr;
        ext->cm_matrix_view->z = nullptr;
        ext->cm_matrix_view->stype = -1; // UPPER in CSR, but LOWER in CSC
        ext->cm_matrix_view->itype = _getCholmodItype<I>();
        ext->cm_matrix_view->xtype = _getCholmodXtype<T>();
        ext->cm_matrix_view->dtype = _getCholmodDtype<T>();
        ext->cm_matrix_view->sorted = 1;
        ext->cm_matrix_view->packed = 1;
    }

    ext->cm_matrix_view->p = a.rows;
    ext->cm_matrix_view->i = a.cols;
    ext->cm_matrix_view->x = a.vals;

    if(ext->stage == 1) ext->stage = 2;
    if(ext->stage == 4) ext->stage = 3;
}

template <typename T, typename I>
void DirectSparseSolver<T, I>::symbolicFactorization(int fixedSuffix)
{
    if(fixedSuffix != 0) eslog::error("symbolicFactorization: dont know what to do with that. TODO\n");

    if(ext->stage != 2) throw std::runtime_error("symbolicFactorization: invalid order of operations in spsolver\n");

    ext->cm_factor_super = _analyze<I>(ext->cm_matrix_view, ext->cm_common);
    ext->factor_nnz = static_cast<I>(ext->cm_common.lnz);

    if(!ext->getfactor_preprocess_lazy) getfactor_preprocess(ext);
    
    ext->stage = 3;
}

template <typename T, typename I>
void DirectSparseSolver<T, I>::numericalFactorization()
{
    if(ext->stage < 3) eslog::error("numericalFactorization: invalid order of operations in spsolver\n");

    _factorize<I>(ext->cm_factor_super, ext->cm_matrix_view, ext->cm_common);

    ext->stage = 4;
}

template <typename T, typename I>
static void DSSsolve(std::unique_ptr<Solver_External_Representation<T,I>> &ext, Vector_Dense<T, I> &rhs, Vector_Dense<T, I> &solution, int sys)
{
    if(ext->stage < 4) eslog::error("solve: invalid order of operations in spsolver\n");

    cholmod_dense cm_rhs;
    cm_rhs.nrow = rhs.size;
    cm_rhs.ncol = 1;
    cm_rhs.d = rhs.size;
    cm_rhs.nzmax = rhs.size;
    cm_rhs.x = rhs.vals;
    cm_rhs.xtype = _getCholmodXtype<T>();
    cm_rhs.dtype = _getCholmodDtype<T>();

    cholmod_dense * cm_sol = _solve<I>(sys, ext->cm_factor_super, &cm_rhs, ext->cm_common);

    solution.resize(cm_sol->nrow);
    std::copy_n(reinterpret_cast<T*>(cm_sol->x), cm_sol->nrow, solution.vals);

    _free<I>(cm_sol, ext->cm_common);
}

template <typename T, typename I>
void DirectSparseSolver<T, I>::solve(Vector_Dense<T, I> &rhs, Vector_Dense<T, I> &solution, int sparsity)
{
    DSSsolve(ext, rhs, solution, CHOLMOD_A);
}

template <typename T, typename I>
void DirectSparseSolver<T, I>::solveForward (Vector_Dense<T, I> &rhs, Vector_Dense<T, I> &solution, int sparsity)
{
    Vector_Dense<T, I> tmp;
    DSSsolve(ext, rhs, solution, CHOLMOD_P);
    DSSsolve(ext, solution, tmp, CHOLMOD_L);
    DSSsolve(ext, tmp, solution, CHOLMOD_P);
}

template <typename T, typename I>
void DirectSparseSolver<T, I>::solveBackward(Vector_Dense<T, I> &rhs, Vector_Dense<T, I> &solution, int sparsity)
{
    Vector_Dense<T, I> tmp;
    DSSsolve(ext, rhs, solution, CHOLMOD_Pt);
    DSSsolve(ext, solution, tmp, CHOLMOD_Lt);
    DSSsolve(ext, tmp, solution, CHOLMOD_Pt);
}

template <typename T, typename I>
void DirectSparseSolver<T, I>::solveDiagonal(Vector_Dense<T, I> &rhs, Vector_Dense<T, I> &solution, int sparsity)
{
    eslog::error("implement solveDiagonal for SuiteSparse.\n");
}

template <typename T, typename I>
static void DSSsolve(std::unique_ptr<Solver_External_Representation<T,I>> &ext, Matrix_Dense<T, I> &rhs, Matrix_Dense<T, I> &solution, int sys)
{
    if(ext->stage < 4) eslog::error("solve: invalid order of operations in spsolver\n");

    cholmod_dense cm_rhs;
    cm_rhs.nrow = rhs.ncols;
    cm_rhs.ncol = rhs.nrows;
    cm_rhs.d = rhs.get_ld();
    cm_rhs.nzmax = cm_rhs.d * rhs.nrows;
    cm_rhs.x = rhs.vals;
    cm_rhs.xtype = _getCholmodXtype<T>();
    cm_rhs.dtype = _getCholmodDtype<T>();

    cholmod_dense * cm_sol = _solve<I>(sys, ext->cm_factor_super, &cm_rhs, ext->cm_common);

    solution.resize(cm_sol->ncol, cm_sol->d);
    solution.ncols = cm_sol->nrow;
    std::copy_n(reinterpret_cast<T*>(cm_sol->x), cm_sol->nzmax, solution.vals);

    _free<I>(cm_sol, ext->cm_common);
}

template <typename T, typename I>
void DirectSparseSolver<T, I>::solve(Matrix_Dense<T, I> &rhs, Matrix_Dense<T, I> &solution, int sparsity)
{
    DSSsolve(ext, rhs, solution, CHOLMOD_A);
}

template <typename T, typename I>
void DirectSparseSolver<T, I>::solveForward (Matrix_Dense<T, I> &rhs, Matrix_Dense<T, I> &solution, int sparsity)
{
    // U
    Matrix_Dense<T, I> tmp;
    DSSsolve(ext, rhs, solution, CHOLMOD_P);
    DSSsolve(ext, solution, tmp, CHOLMOD_L);
    DSSsolve(ext, tmp, solution, CHOLMOD_P);
}

template <typename T, typename I>
void DirectSparseSolver<T, I>::solveBackward(Matrix_Dense<T, I> &rhs, Matrix_Dense<T, I> &solution, int sparsity)
{
    // L
    Matrix_Dense<T, I> tmp;
    DSSsolve(ext, rhs, solution, CHOLMOD_Pt);
    DSSsolve(ext, solution, tmp, CHOLMOD_Lt);
    DSSsolve(ext, tmp, solution, CHOLMOD_Pt);
}

template <typename T, typename I>
void DirectSparseSolver<T, I>::solveDiagonal(Matrix_Dense<T, I> &rhs, Matrix_Dense<T, I> &solution, int sparsity)
{
    eslog::error("implement solveDiagonal for SuiteSparse.\n");
}

template <typename T, typename I>
I DirectSparseSolver<T, I>::getMatrixSize()
{
    return ext->cm_matrix_view->nrow;
}

template <typename T, typename I>
I DirectSparseSolver<T, I>::getMatrixNnz()
{
    return ext->cm_matrix_view->nzmax;
}

template <typename T, typename I>
I DirectSparseSolver<T, I>::getFactorNnz()
{
    if(ext->stage < 3) eslog::error("getFactorNnz: invalid order of operations in spsolver\n");

    // https://github.com/DrTimothyAldenDavis/SuiteSparse/issues/523

    return ext->factor_nnz;
}

template <typename T, typename I>
template <typename A>
void DirectSparseSolver<T, I>::getFactorL(Matrix_CSR<T,I,A> &/*L*/, bool /*copyPattern*/, bool /*copyValues*/)
{
    eslog::error("L factor is not provided\n");
}

template <typename T, typename I>
template <typename A>
void DirectSparseSolver<T, I>::getFactorU(Matrix_CSR<T,I,A> &U, bool copyPattern, bool copyValues)
{
    if(ext->stage < 3) eslog::error("getFactorU: invalid order of operations in spsolver\n");
    if(copyValues && ext->stage < 4) eslog::error("getFactorU: invalid order of operations in spsolver\n");
    if((size_t)U.nrows != ext->cm_factor_super->n || (size_t)U.ncols != ext->cm_factor_super->n || U.nnz != getFactorNnz()) eslog::error("getFactorU: output matrix has wrong dimensions\n");

    if(!ext->getfactor_preprocess_done) getfactor_preprocess(ext);

    if(copyPattern) {
        std::copy_n(static_cast<I*>(ext->cm_factor_simpl->p), ext->cm_factor_simpl->n+1, U.rows);
        std::copy_n(static_cast<I*>(ext->cm_factor_simpl->i), ext->cm_factor_simpl->nzmax, U.cols);
    }
    if(copyValues) {
        for(I i = 0; i < ext->map_simpl_super.size; i++) {
            U.vals[i] = reinterpret_cast<T*>(ext->cm_factor_super->x)[ext->map_simpl_super.vals[i]];
        }
    }
}

template <typename T, typename I>
void DirectSparseSolver<T, I>::getPermutation(Permutation<I> &perm)
{
    if(ext->stage < 3) eslog::error("getPermutation: invalid order of operations in spsolver\n");

    perm.resize(ext->cm_factor_super->n);

    std::copy_n(static_cast<I*>(ext->cm_factor_super->Perm), ext->cm_factor_super->n, perm.dst_to_src);

    if(ext->cm_factor_super->IPerm != nullptr) {
        std::copy_n(static_cast<I*>(ext->cm_factor_super->IPerm), ext->cm_factor_super->n, perm.dst_to_src);
    }
    else {
        perm.invert(perm.dst_to_src, perm.src_to_dst);
    }
}

template <typename T, typename I>
void DirectSparseSolver<T, I>::getPermutation(Vector_Dense<I> &perm)
{
    if(ext->stage < 3) eslog::error("getPermutation: invalid order of operations in spsolver\n");

    perm.resize(ext->cm_factor_super->n);

    std::copy_n(static_cast<I*>(ext->cm_factor_super->Perm), ext->cm_factor_super->n, perm.vals);
}

template <typename T, typename I>
void DirectSparseSolver<T, I>::getSC(Matrix_Dense<T,I> &sc)
{
    // todo: can I use the existing factor so i don't have to factorize again?
    // computes the schur complement S = A22 - A21 * A11^{-1} * A12, where A = [A11, A12; A21, A22]

    if(ext->stage < 2) eslog::error("getSC: invalid order of operations in spsolver\n");

    I size_sc = sc.nrows;
    I size = ext->matrix->nrows;
    I size_A11 = size - size_sc;

    Matrix_CSR<T, I> A11_sp;
    Matrix_CSR<T, I> A21t_sp; // = A12c_sp
    Matrix_Dense<T, I> A22t_dn;
    Matrix_Dense<T, I> A12t_dn;
    SpBLAS<Matrix_CSR, T, I>::submatrix(*ext->matrix, A11_sp, 0, size_A11, 0, size_A11);
    SpBLAS<Matrix_CSR, T, I>::submatrix(*ext->matrix, A21t_sp, 0, size_A11, size_A11, size, false, true); // = A12c_sp
    SpBLAS<Matrix_CSR, T, I>::submatrix(*ext->matrix, A22t_dn, size_A11, size, size_A11, size, true, false, true);
    SpBLAS<Matrix_CSR, T, I>::submatrix(*ext->matrix, A12t_dn, 0, size_A11, size_A11, size, true, false, true);

    cholmod_sparse *cm_A11_sp = nullptr;
    cholmod_sparse *cm_A21_sp = nullptr;
    cholmod_dense *cm_A22_dn = nullptr;
    cholmod_dense *cm_A12_dn = nullptr;
    cholmod_factor *cm_L = nullptr;
    cholmod_dense *cm_A11iA12_dn = nullptr;

    setSymmetric(cm_A11_sp, A11_sp);
    updateSymmetric(cm_A11_sp, A11_sp);
    setSymmetric(cm_A21_sp, A21t_sp);
    updateSymmetric(cm_A21_sp, A21t_sp);
    update(cm_A22_dn, A22t_dn);
    update(cm_A12_dn, A12t_dn);

    double alpha[2] = {-1,0};
    double beta[2] = {1,0};

    cm_L = _analyze<I>(cm_A11_sp, ext->cm_common);
    _factorize<I>(cm_L, cm_A11_sp, ext->cm_common);
    cm_A11iA12_dn = _solve<I>(CHOLMOD_A, cm_L, cm_A12_dn, ext->cm_common);
    _apply<I>(cm_A22_dn, cm_A21_sp, cm_A11iA12_dn, alpha, beta, ext->cm_common);

    if constexpr (std::is_same_v<T,double>) { sc.type = Matrix_Type::REAL_SYMMETRIC_POSITIVE_DEFINITE; }
    if constexpr (std::is_same_v<T,std::complex<double>>) { sc.type = Matrix_Type::COMPLEX_HERMITIAN_POSITIVE_DEFINITE; }
    sc.shape = Matrix_Shape::UPPER;
    sc.resize(A22t_dn);
    for(I r = 0, i = 0; r < sc.nrows; ++r) {
        for(I c = r; c < sc.ncols; ++c, ++i) {
            sc.vals[i] = A22t_dn.vals[r * sc.ncols + c];
        }
    }

    delete cm_A11_sp;
    delete cm_A21_sp;
    delete cm_A22_dn;
    delete cm_A12_dn;
    _free<I>(cm_L, ext->cm_common);
    _free<I>(cm_A11iA12_dn, ext->cm_common);
}

}

#include "math/wrappers/math.spsolver.inst.hpp"

#endif

