
#include "math/math.h"
#include "esinfo/eslog.h"
#include "gpu/gpu_management.h"
#include "basis/utilities/utils.h"

#include <complex>

#ifdef HAVE_SUITESPARSE
#ifdef ESPRESO_USE_WRAPPER_SPSOLVER_SUITESPARSE

#include "math/wrappers/math.spsolver.h"
#include "w.suitesparse.cholmod.h"
#include "w.suitesparse.umfpack.h"

namespace espreso {

template<typename T, typename I>
struct Cholmod_Solver_External_Representation {
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
struct Umfpack_Solver_External_Representation {
    const Matrix_CSR<T, I> * matrix = nullptr;
    Matrix_CSC<T, I> full;

    void *symbolic, *numeric;
    double control[UMFPACK_CONTROL], info[UMFPACK_INFO];
};

template<typename T, typename I>
struct Solver_External_Representation {
    enum struct SOLVER {
        CHOLMOD,
        UMFPACK
    };

    SOLVER solver;
    Cholmod_Solver_External_Representation<T, I> cholmod;
    Umfpack_Solver_External_Representation<T, I> umfpack;
};

static void check(double info[])
{
    switch ((int)info[UMFPACK_STATUS]) {
    case UMFPACK_OK: break;
    case UMFPACK_ERROR_n_nonpositive:    eslog::error("UMFPACK_ERROR: n_nonpositive.\n"); break;
    case UMFPACK_ERROR_invalid_matrix:   eslog::error("UMFPACK_ERROR: invalid_matrix.\n"); break;
    case UMFPACK_ERROR_out_of_memory:    eslog::error("UMFPACK_ERROR: out_of_memory.\n"); break;
    case UMFPACK_ERROR_argument_missing: eslog::error("UMFPACK_ERROR: argument_missing.\n"); break;
    case UMFPACK_ERROR_internal_error:   eslog::error("UMFPACK_ERROR: internal_error.\n"); break;
    }
}

template<typename T, typename I>
void getfactor_preprocess(std::unique_ptr<Solver_External_Representation<T,I>> & ext)
{
    switch (ext->solver) {
    case Solver_External_Representation<T, I>::SOLVER::CHOLMOD: {
        if (ext->cholmod.cm_factor_super->xsize > utils::get_max_val_no_precision_loss_in_fp<T>()) eslog::error("symbolicFactorization: factor nnz too large for my super->simpl map\n");

        ext->cholmod.cm_factor_simpl = _copyFactor<I>(ext->cholmod.cm_factor_super, ext->cholmod.cm_common);

        _changeFactor<I>(_getCholmodXtype<T>(), true, true, true, true, ext->cholmod.cm_factor_simpl, ext->cholmod.cm_common);

        for(size_t i = 0; i < ext->cholmod.cm_factor_simpl->xsize; i++) {
            reinterpret_cast<T*>(ext->cholmod.cm_factor_simpl->x)[i] = static_cast<T>(i);
        }

        _changeFactor<I>(_getCholmodXtype<T>(), true, false, true, true, ext->cholmod.cm_factor_simpl, ext->cholmod.cm_common);

        if (ext->cholmod.zerodrop == 'D') {
            _resymbol<I>(ext->cholmod.cm_matrix_view, nullptr, 0, 1, ext->cholmod.cm_factor_simpl, ext->cholmod.cm_common);
        }

        if (ext->cholmod.cm_factor_simpl->nzmax != (size_t)ext->cholmod.factor_nnz) {
            eslog::error("getfactor_preprocess: some weird error in cholmod\n");
        }

        ext->cholmod.map_simpl_super.resize(ext->cholmod.factor_nnz);
        for(I i = 0; i < ext->cholmod.map_simpl_super.size; i++) {
            ext->cholmod.map_simpl_super.vals[i] = static_cast<I>(std::real(reinterpret_cast<T*>(ext->cholmod.cm_factor_simpl->x)[i]));
        }

        ext->cholmod.getfactor_preprocess_done = true;
    } break;
    case Solver_External_Representation<T, I>::SOLVER::UMFPACK: {

    } break;
    }
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
    return false;
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

    ext->cholmod.zerodrop = 'D'; // Drop, Keep
    ext->cholmod.getfactor_preprocess_lazy = true;

    _start<I>(ext->cholmod.cm_common);
    ext->cholmod.cm_common.final_ll = 1;
    ext->cholmod.cm_common.nthreads_max = 1;
    ext->cholmod.cm_common.nmethods = 1;
    ext->cholmod.cm_common.method[0].ordering = CHOLMOD_METIS;
    ext->cholmod.cm_common.itype = _getCholmodItype<I>();
    ext->cholmod.cm_common.supernodal = CHOLMOD_SUPERNODAL;

    ext->cholmod.stage = 1;

    ext->umfpack.matrix = nullptr;
    ext->umfpack.symbolic = nullptr;
    ext->umfpack.numeric = nullptr;
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
    switch (ext->solver) {
    case Solver_External_Representation<T, I>::SOLVER::CHOLMOD: {
        delete ext->cholmod.cm_matrix_view;

        if (ext->cholmod.cm_factor_super != nullptr) _free<I>(ext->cholmod.cm_factor_super, ext->cholmod.cm_common);

        if (ext->cholmod.cm_factor_simpl != nullptr) _free<I>(ext->cholmod.cm_factor_simpl, ext->cholmod.cm_common);

        _finish<I>(ext->cholmod.cm_common);
    } break;
    case Solver_External_Representation<T, I>::SOLVER::UMFPACK: {
        if (ext->umfpack.symbolic != nullptr) umfpack_di_free_symbolic(&ext->umfpack.symbolic);
        if (ext->umfpack.numeric != nullptr) umfpack_di_free_numeric(&ext->umfpack.numeric);
    }
    }
}

template <typename T, typename I>
void DirectSparseSolver<T, I>::commit(const Matrix_CSR<T,I> &a)
{
    if (a.nrows != a.ncols) eslog::error("commit: matrix has to be square\n");
    if (a.rows[0] != 0) eslog::error("commit: CSR matrix has to use 0-based indexing\n");

    if (a.type == Matrix_Type::REAL_SYMMETRIC_POSITIVE_DEFINITE || a.type == Matrix_Type::COMPLEX_HERMITIAN_POSITIVE_DEFINITE) {
        ext->solver = Solver_External_Representation<T, I>::SOLVER::CHOLMOD;
    } else {
        ext->solver = Solver_External_Representation<T, I>::SOLVER::UMFPACK;
    }

    switch (ext->solver) {
    case Solver_External_Representation<T, I>::SOLVER::CHOLMOD: {
        if (ext->cholmod.stage < 1) eslog::error("commit: invalid order of operations in spsolver\n");

        ext->cholmod.matrix = &a;
        if (ext->cholmod.cm_matrix_view == nullptr) {
            ext->cholmod.cm_matrix_view = new cholmod_sparse();
            ext->cholmod.cm_matrix_view->nrow = a.ncols;
            ext->cholmod.cm_matrix_view->ncol = a.nrows;
            ext->cholmod.cm_matrix_view->nzmax = a.nnz;
            ext->cholmod.cm_matrix_view->nz = nullptr;
            ext->cholmod.cm_matrix_view->z = nullptr;
            ext->cholmod.cm_matrix_view->stype = -1; // UPPER in CSR, but LOWER in CSC
            ext->cholmod.cm_matrix_view->itype = _getCholmodItype<I>();
            ext->cholmod.cm_matrix_view->xtype = _getCholmodXtype<T>();
            ext->cholmod.cm_matrix_view->dtype = _getCholmodDtype<T>();
            ext->cholmod.cm_matrix_view->sorted = 1;
            ext->cholmod.cm_matrix_view->packed = 1;
        }

        ext->cholmod.cm_matrix_view->p = a.rows;
        ext->cholmod.cm_matrix_view->i = a.cols;
        ext->cholmod.cm_matrix_view->x = a.vals;

        if (ext->cholmod.stage == 1) ext->cholmod.stage = 2;
        if (ext->cholmod.stage == 4) ext->cholmod.stage = 3;

    } break;
    case Solver_External_Representation<T, I>::SOLVER::UMFPACK: {
        ext->umfpack.matrix = &a;
        _defaults<I>(ext->umfpack.control);
        switch (a.type) {
        case Matrix_Type::UNSET_INVALID_NONE:
            eslog::error("Invalid/unset matrix type\n");
        case Matrix_Type::REAL_SYMMETRIC_INDEFINITE:
        case Matrix_Type::REAL_SYMMETRIC_POSITIVE_DEFINITE:
        case Matrix_Type::COMPLEX_SYMMETRIC:
        case Matrix_Type::COMPLEX_HERMITIAN_INDEFINITE:
        case Matrix_Type::COMPLEX_HERMITIAN_POSITIVE_DEFINITE:
            ext->umfpack.control[UMFPACK_STRATEGY] = UMFPACK_STRATEGY_SYMMETRIC;
            break;
        case Matrix_Type::REAL_NONSYMMETRIC:
        case Matrix_Type::REAL_STRUCTURALLY_SYMMETRIC:
        case Matrix_Type::COMPLEX_NONSYMMETRIC:
        case Matrix_Type::COMPLEX_STRUCTURALLY_SYMMETRIC:
            ext->umfpack.control[UMFPACK_STRATEGY] = UMFPACK_STRATEGY_UNSYMMETRIC;
            break;
        }

        if (a.nrows == 0 || a.ncols == 0) {
            return;
        }

        switch (a.shape) { // UMFPACK need fully stored matrix in CSC
        case Matrix_Shape::LOWER: eslog::error("implement DirectSparseSolver for LOWER triangular matrices.\n"); break;
        case Matrix_Shape::FULL:  ext->umfpack.full.pattern(a); break;
        case Matrix_Shape::UPPER: {
            Matrix_CSC<T, I> &full = ext->umfpack.full;
            full.resize(a.nrows, a.ncols, 2 * (a.nnz - a.nrows) + a.nrows);
            std::fill(full.cols, full.cols + full.ncols, 0);
            for (I r = 0; r < a.nrows; ++r) {
                full.cols[r] += a.rows[r + 1] - a.rows[r];
                for (I c = a.rows[r] + 1; c < a.rows[r + 1]; ++c) {
                    ++full.cols[a.cols[c - Indexing::CSR]];
                }
            }
            utils::sizesToOffsets(full.cols, full.cols + full.ncols + 1);
            for (I r = 0; r < a.nrows; ++r) {
                for (I c = a.rows[r]; c < a.rows[r + 1]; ++c) {
                    // copy whole column
                    full.rows[full.cols[r]++] = a.cols[c - Indexing::CSR];
                    // copy prefix of upper triangle
                    if (c != a.rows[r]) { // skip the diagonal
                        full.rows[full.cols[a.cols[c - Indexing::CSR]]++] = r;
                    }
                }
            }
            for (int r = 0, next = 0; r <= full.ncols; ++r) { // reset row pointers
                int current = next;
                next = full.cols[r];
                full.cols[r] = current;
            }
        } break;
        }
    } break;
    }
}

template <typename T, typename I>
void DirectSparseSolver<T, I>::symbolicFactorization()
{
    switch (ext->solver) {
    case Solver_External_Representation<T, I>::SOLVER::CHOLMOD: {
        if (ext->cholmod.cm_matrix_view->nrow == 0 || ext->cholmod.cm_matrix_view->ncol == 0) return;
        if (ext->cholmod.stage != 2) eslog::error("symbolicFactorization: invalid order of operations in spsolver\n");

        ext->cholmod.cm_factor_super = _analyze<I>(ext->cholmod.cm_matrix_view, ext->cholmod.cm_common);

        if(ext->cholmod.cm_common.lnz > (double)std::numeric_limits<I>::max()) eslog::error("symbolicFactorization: factor nnz too large for the used integer type\n");
        ext->cholmod.factor_nnz = static_cast<I>(ext->cholmod.cm_common.lnz);

        if (!ext->cholmod.getfactor_preprocess_lazy) getfactor_preprocess(ext);

        ext->cholmod.stage = 3;
    } break;

    case Solver_External_Representation<T, I>::SOLVER::UMFPACK: {
        if (ext->umfpack.full.nrows == 0 || ext->umfpack.full.ncols == 0) return;
        if (ext->umfpack.symbolic != nullptr) eslog::error("symbolicFactorization: invalid order of operations in spsolver\n");

        _symbolic<T, I>(ext->umfpack.full, &ext->umfpack.symbolic, ext->umfpack.control, ext->umfpack.info);
        check(ext->umfpack.info);
    } break;
    }
}

template <typename T, typename I>
void DirectSparseSolver<T, I>::numericalFactorization()
{
    switch (ext->solver) {
    case Solver_External_Representation<T, I>::SOLVER::CHOLMOD: {
        if (ext->cholmod.cm_matrix_view->nrow == 0 || ext->cholmod.cm_matrix_view->ncol == 0) return;
        if (ext->cholmod.stage < 3) eslog::error("numericalFactorization: invalid order of operations in spsolver\n");

        _factorize<I>(ext->cholmod.cm_factor_super, ext->cholmod.cm_matrix_view, ext->cholmod.cm_common);

        ext->cholmod.stage = 4;
    } break;

    case Solver_External_Representation<T, I>::SOLVER::UMFPACK: {
        if (ext->umfpack.full.nrows == 0 || ext->umfpack.full.ncols == 0) return;
        if (ext->umfpack.symbolic == nullptr) eslog::error("numericalFactorization: invalid order of operations in spsolver\n");

        switch (ext->umfpack.matrix->shape) { // UMFPACK need fully stored matrix in CSC
        case Matrix_Shape::LOWER: eslog::error("implement DirectSparseSolver for LOWER triangular matrices.\n"); break;
        case Matrix_Shape::FULL: {
            eslog::error("check CORRECTNESS.\n");
            const Matrix_CSR<T, I> *a = ext->umfpack.matrix;
            Matrix_CSC<T, I> &full = ext->umfpack.full;

            for (I r = 0; r < a->nrows; ++r) {
                for (I c = a->rows[r]; c < a->rows[r + 1]; ++c) {
                    full.vals[full.cols[a->cols[c - Indexing::CSR]]++] = a->vals[c - Indexing::CSR];
                }
            }
            for (int r = 0, next = 0; r <= full.ncols; ++r) { // reset row pointers
                int current = next;
                next = full.cols[r];
                full.cols[r] = current;
            }
        } break;
        case Matrix_Shape::UPPER: {
            const Matrix_CSR<T, I> *a = ext->umfpack.matrix;
            Matrix_CSC<T, I> &full = ext->umfpack.full;

            for (I r = 0; r < a->nrows; ++r) {
                for (I c = a->rows[r]; c < a->rows[r + 1]; ++c) {
                    full.vals[full.cols[r]++] = a->vals[c - Indexing::CSR];
                    if (c != a->rows[r]) { // skip the diagonal
                        full.vals[full.cols[a->cols[c - Indexing::CSR]]++] = a->vals[c - Indexing::CSR];
                    }
                }
            }
            for (int r = 0, next = 0; r <= full.ncols; ++r) { // reset row pointers
                int current = next;
                next = full.cols[r];
                full.cols[r] = current;
            }
        } break;
        }

        _numeric<T, I>(ext->umfpack.full, ext->umfpack.symbolic, &ext->umfpack.numeric, ext->umfpack.control, ext->umfpack.info);
        check(ext->umfpack.info);
    } break;
    }
}

template <typename T, typename I>
static void DSSsolve(std::unique_ptr<Solver_External_Representation<T,I>> &ext, const Vector_Dense<T, I> &rhs, Vector_Dense<T, I> &solution, int sys)
{
    if (rhs.size == 0 || solution.size == 0) return;
    if(rhs.size != ext->cholmod.matrix->nrows) eslog::error("wrong rhs size\n");
    if(solution.size != ext->cholmod.matrix->nrows) eslog::error("wrong solution size\n");

    switch (ext->solver) {
    case Solver_External_Representation<T, I>::SOLVER::CHOLMOD: {
        if (ext->cholmod.cm_matrix_view->nrow == 0 || ext->cholmod.cm_matrix_view->ncol == 0) return;
        if (ext->cholmod.stage < 4) eslog::error("solve: invalid order of operations in spsolver\n");

        cholmod_dense cm_rhs;
        cm_rhs.nrow = rhs.size;
        cm_rhs.ncol = 1;
        cm_rhs.d = rhs.size;
        cm_rhs.nzmax = rhs.size;
        cm_rhs.x = rhs.vals;
        cm_rhs.xtype = _getCholmodXtype<T>();
        cm_rhs.dtype = _getCholmodDtype<T>();

        cholmod_dense * cm_sol = _solve<I>(sys, ext->cholmod.cm_factor_super, &cm_rhs, ext->cholmod.cm_common);

        std::copy_n(reinterpret_cast<T*>(cm_sol->x), cm_sol->nrow, solution.vals);

        _free<I>(cm_sol, ext->cholmod.cm_common);
    } break;
    case Solver_External_Representation<T, I>::SOLVER::UMFPACK: {
        if (ext->umfpack.full.nrows == 0 || ext->umfpack.full.ncols == 0) return;
        if (ext->umfpack.numeric == nullptr) eslog::error("solve: invalid order of operations in spsolver\n");

        _solve<T, I>(UMFPACK_A, ext->umfpack.full, solution.vals, rhs.vals, ext->umfpack.numeric, ext->umfpack.control, ext->umfpack.info);
        check(ext->umfpack.info);
    } break;
    }
}

template <typename T, typename I>
void DirectSparseSolver<T, I>::solve(const Vector_Dense<T, I> &rhs, Vector_Dense<T, I> &solution)
{
    int sys = 0;
    switch (ext->solver) {
    case Solver_External_Representation<T, I>::SOLVER::CHOLMOD: sys = CHOLMOD_A; break;
    case Solver_External_Representation<T, I>::SOLVER::UMFPACK: sys = UMFPACK_A; break;
    }
    DSSsolve(ext, rhs, solution, sys);
}

template <typename T, typename I>
void DirectSparseSolver<T, I>::solveForward (const Vector_Dense<T, I> &rhs, Vector_Dense<T, I> &solution)
{
    Vector_Dense<T, I> tmp;
    switch (ext->solver) {
    case Solver_External_Representation<T, I>::SOLVER::CHOLMOD:
        DSSsolve(ext, rhs, solution, CHOLMOD_P);
        DSSsolve(ext, solution, tmp, CHOLMOD_L);
        DSSsolve(ext, tmp, solution, CHOLMOD_P);
    break;
    case Solver_External_Representation<T, I>::SOLVER::UMFPACK:
        eslog::error("implement me\n");
    break;
    }
}

template <typename T, typename I>
void DirectSparseSolver<T, I>::solveBackward(const Vector_Dense<T, I> &rhs, Vector_Dense<T, I> &solution)
{
    Vector_Dense<T, I> tmp;
    switch (ext->solver) {
    case Solver_External_Representation<T, I>::SOLVER::CHOLMOD:
        DSSsolve(ext, rhs, solution, CHOLMOD_Pt);
        DSSsolve(ext, solution, tmp, CHOLMOD_Lt);
        DSSsolve(ext, tmp, solution, CHOLMOD_Pt);
    break;
    case Solver_External_Representation<T, I>::SOLVER::UMFPACK:
        eslog::error("implement me\n");
    break;
    }
}

template <typename T, typename I>
void DirectSparseSolver<T, I>::solveDiagonal(const Vector_Dense<T, I> &rhs, Vector_Dense<T, I> &solution)
{
    eslog::error("implement solveDiagonal for SuiteSparse.\n");
}

template <typename T, typename I>
static void DSSsolve(std::unique_ptr<Solver_External_Representation<T,I>> &ext, const Matrix_Dense<T, I> &rhs, Matrix_Dense<T, I> &solution, int sys)
{
    if (rhs.nrows == 0 || rhs.ncols == 0) return;

    switch (ext->solver) {
    case Solver_External_Representation<T, I>::SOLVER::CHOLMOD: {
        if (ext->cholmod.cm_matrix_view->nrow == 0 || ext->cholmod.cm_matrix_view->ncol == 0) return;
        if (ext->cholmod.stage < 4) eslog::error("solve: invalid order of operations in spsolver\n");

        cholmod_dense cm_rhs;
        cm_rhs.nrow = rhs.ncols;
        cm_rhs.ncol = rhs.nrows;
        cm_rhs.d = rhs.get_ld();
        cm_rhs.nzmax = cm_rhs.d * rhs.nrows;
        cm_rhs.x = rhs.vals;
        cm_rhs.xtype = _getCholmodXtype<T>();
        cm_rhs.dtype = _getCholmodDtype<T>();

        cholmod_dense * cm_sol = _solve<I>(sys, ext->cholmod.cm_factor_super, &cm_rhs, ext->cholmod.cm_common);

        solution.resize(cm_sol->ncol, cm_sol->d);
        solution.ncols = cm_sol->nrow;
        std::copy_n(reinterpret_cast<T*>(cm_sol->x), cm_sol->nzmax, solution.vals);

        _free<I>(cm_sol, ext->cholmod.cm_common);
    } break;
    case Solver_External_Representation<T, I>::SOLVER::UMFPACK: {
        if (ext->umfpack.full.nrows == 0 || ext->umfpack.full.ncols == 0) return;
        if (ext->umfpack.numeric == nullptr) eslog::error("solve: invalid order of operations in spsolver\n");

        solution.resize(rhs.nrows, rhs.ncols);
        for (int c = 0; c < rhs.nrows; ++c) {
            _solve<T, I>(sys, ext->umfpack.full, solution.vals + c * solution.ncols, rhs.vals + c * rhs.ncols, ext->umfpack.numeric, ext->umfpack.control, ext->umfpack.info);
            check(ext->umfpack.info);
        }
    } break;
    }
}

template <typename T, typename I>
void DirectSparseSolver<T, I>::solve(const Matrix_Dense<T, I> &rhs, Matrix_Dense<T, I> &solution)
{
    int sys = 0;
    switch (ext->solver) {
    case Solver_External_Representation<T, I>::SOLVER::CHOLMOD: sys = CHOLMOD_A; break;
    case Solver_External_Representation<T, I>::SOLVER::UMFPACK: sys = UMFPACK_A; break;
    }
    DSSsolve(ext, rhs, solution, sys);
}

template <typename T, typename I>
void DirectSparseSolver<T, I>::solveForward (const Matrix_Dense<T, I> &rhs, Matrix_Dense<T, I> &solution)
{
    // U
    Matrix_Dense<T, I> tmp;
    switch (ext->solver) {
    case Solver_External_Representation<T, I>::SOLVER::CHOLMOD:
        DSSsolve(ext, rhs, solution, CHOLMOD_P);
        DSSsolve(ext, solution, tmp, CHOLMOD_L);
        DSSsolve(ext, tmp, solution, CHOLMOD_P);
    break;
    case Solver_External_Representation<T, I>::SOLVER::UMFPACK:
        eslog::error("implement me\n");
    break;
    }
}

template <typename T, typename I>
void DirectSparseSolver<T, I>::solveBackward(const Matrix_Dense<T, I> &rhs, Matrix_Dense<T, I> &solution)
{
    // L
    Matrix_Dense<T, I> tmp;
    switch (ext->solver) {
    case Solver_External_Representation<T, I>::SOLVER::CHOLMOD:
        DSSsolve(ext, rhs, solution, CHOLMOD_Pt);
        DSSsolve(ext, solution, tmp, CHOLMOD_Lt);
        DSSsolve(ext, tmp, solution, CHOLMOD_Pt);
    break;
    case Solver_External_Representation<T, I>::SOLVER::UMFPACK:
        eslog::error("implement me\n");
    break;
    }
}

template <typename T, typename I>
void DirectSparseSolver<T, I>::solveDiagonal(const Matrix_Dense<T, I> &rhs, Matrix_Dense<T, I> &solution)
{
    eslog::error("implement solveDiagonal for SuiteSparse.\n");
}

template <typename T, typename I>
I DirectSparseSolver<T, I>::getMatrixSize()
{
    switch (ext->solver) {
    case Solver_External_Representation<T, I>::SOLVER::CHOLMOD: return ext->cholmod.cm_matrix_view->nrow;
    case Solver_External_Representation<T, I>::SOLVER::UMFPACK: return ext->umfpack.matrix->nrows;
    }
    return 0;
}

template <typename T, typename I>
I DirectSparseSolver<T, I>::getMatrixNnz()
{
    switch (ext->solver) {
    case Solver_External_Representation<T, I>::SOLVER::CHOLMOD: return ext->cholmod.cm_matrix_view->nzmax;
    case Solver_External_Representation<T, I>::SOLVER::UMFPACK: return ext->umfpack.matrix->nnz;
    }
    return 0;
}

template <typename T, typename I>
I DirectSparseSolver<T, I>::getFactorNnz()
{
    switch (ext->solver) {
    case Solver_External_Representation<T, I>::SOLVER::CHOLMOD: {
        if (ext->cholmod.stage < 3) eslog::error("getFactorNnz: invalid order of operations in spsolver\n");

        // https://github.com/DrTimothyAldenDavis/SuiteSparse/issues/523

        return ext->cholmod.factor_nnz;
    }

    case Solver_External_Representation<T, I>::SOLVER::UMFPACK: {
        return ext->umfpack.info[UMFPACK_LNZ_ESTIMATE] + ext->umfpack.info[UMFPACK_UNZ_ESTIMATE];
    }
    };
    return 0;
}

template <typename T, typename I>
template <typename A>
void DirectSparseSolver<T, I>::getFactorL(Matrix_CSR<T,I,A> &/*L*/, bool /*copyPattern*/, bool /*copyValues*/)
{
    switch (ext->solver) {
    case Solver_External_Representation<T, I>::SOLVER::CHOLMOD: {
        eslog::error("L factor is not provided\n");
    } break;
    case Solver_External_Representation<T, I>::SOLVER::UMFPACK: {
        eslog::error("implement DirectSparseSolver<T, I>::getFactorL\n");
    } break;
    }
}

template <typename T, typename I>
template <typename A>
void DirectSparseSolver<T, I>::getFactorU(Matrix_CSR<T,I,A> &U, bool copyPattern, bool copyValues)
{
    switch (ext->solver) {
    case Solver_External_Representation<T, I>::SOLVER::CHOLMOD: {
        if (ext->cholmod.stage < 3) eslog::error("getFactorU: invalid order of operations in spsolver\n");
        if (copyValues && ext->cholmod.stage < 4) eslog::error("getFactorU: invalid order of operations in spsolver\n");
        if ((size_t)U.nrows != ext->cholmod.cm_factor_super->n || (size_t)U.ncols != ext->cholmod.cm_factor_super->n || U.nnz != getFactorNnz()) eslog::error("getFactorU: output matrix has wrong dimensions\n");

        if (!ext->cholmod.getfactor_preprocess_done) getfactor_preprocess(ext);

        if (copyPattern) {
            std::copy_n(static_cast<I*>(ext->cholmod.cm_factor_simpl->p), ext->cholmod.cm_factor_simpl->n+1, U.rows);
            std::copy_n(static_cast<I*>(ext->cholmod.cm_factor_simpl->i), ext->cholmod.cm_factor_simpl->nzmax, U.cols);
        }
        if (copyValues) {
            for(I i = 0; i < ext->cholmod.map_simpl_super.size; i++) {
                U.vals[i] = reinterpret_cast<T*>(ext->cholmod.cm_factor_super->x)[ext->cholmod.map_simpl_super.vals[i]];
            }
        }
    } break;

    case Solver_External_Representation<T, I>::SOLVER::UMFPACK: {
        eslog::error("implement DirectSparseSolver<T, I>::getFactorU\n");
    }
    }
}

template <typename T, typename I>
void DirectSparseSolver<T, I>::getPermutation(Permutation<I> &perm)
{
    switch (ext->solver) {
    case Solver_External_Representation<T, I>::SOLVER::CHOLMOD: {
        if (ext->cholmod.stage < 3) eslog::error("getPermutation: invalid order of operations in spsolver\n");

        if((size_t)perm.size != ext->cholmod.cm_factor_super->n) {
            perm.resize(ext->cholmod.cm_factor_super->n);
        }

        std::copy_n(static_cast<I*>(ext->cholmod.cm_factor_super->Perm), ext->cholmod.cm_factor_super->n, perm.dst_to_src);

        if (ext->cholmod.cm_factor_super->IPerm != nullptr) {
            std::copy_n(static_cast<I*>(ext->cholmod.cm_factor_super->IPerm), ext->cholmod.cm_factor_super->n, perm.dst_to_src);
        }
        else {
            perm.invert(perm.dst_to_src, perm.src_to_dst);
        }
    } break;
    case Solver_External_Representation<T, I>::SOLVER::UMFPACK: {
        eslog::error("implement DirectSparseSolver<T, I>::getPermutation\n");
    } break;
    }
}

template <typename T, typename I>
void DirectSparseSolver<T, I>::getPermutation(Vector_Dense<I> &perm)
{
    switch (ext->solver) {
    case Solver_External_Representation<T, I>::SOLVER::CHOLMOD: {
        if (ext->cholmod.stage < 3) eslog::error("getPermutation: invalid order of operations in spsolver\n");

        if((size_t)perm.size != ext->cholmod.cm_factor_super->n) {
            perm.resize(ext->cholmod.cm_factor_super->n);
        }

        std::copy_n(static_cast<I*>(ext->cholmod.cm_factor_super->Perm), ext->cholmod.cm_factor_super->n, perm.vals);
    } break;
    case Solver_External_Representation<T, I>::SOLVER::UMFPACK: {
        eslog::error("implement DirectSparseSolver<T, I>::getPermutation\n");
    } break;
    }
}

template <typename T, typename I>
void DirectSparseSolver<T, I>::getSC(Matrix_Dense<T,I> &sc, std::vector<I> &indices, bool symmetric_packed)
{
    // todo: can I use the existing factor so i don't have to factorize again?
    // computes the schur complement S = A22 - A21 * A11^{-1} * A12, where A = [A11, A12; A21, A22]

    eslog::error("suitesparse spsolver getSC: not implemented\n");

    // generalized and moved to Dirichlet preconditioner
//    switch (ext->solver) {
//    case Solver_External_Representation<T, I>::SOLVER::CHOLMOD: {
//        if (ext->cholmod.stage < 2) eslog::error("getSC: invalid order of operations in spsolver\n");
//
//        I size_sc = sc.nrows;
//        I size = ext->cholmod.matrix->nrows;
//        I size_A11 = size - size_sc;
//
//        Matrix_CSR<T, I> A11_sp;
//        Matrix_CSR<T, I> A21t_sp; // = A12c_sp
//        Matrix_Dense<T, I> A22t_dn;
//        Matrix_Dense<T, I> A12t_dn;
//        math::spblas::submatrix(*ext->cholmod.matrix, A11_sp, 0, size_A11, 0, size_A11);
//        math::spblas::submatrix(*ext->cholmod.matrix, A21t_sp, 0, size_A11, size_A11, size, false, true); // = A12c_sp
//        math::spblas::submatrix(*ext->cholmod.matrix, A22t_dn, size_A11, size, size_A11, size, true, false, true);
//        math::spblas::submatrix(*ext->cholmod.matrix, A12t_dn, 0, size_A11, size_A11, size, true, false, true);
//
//        cholmod_sparse *cm_A11_sp = nullptr;
//        cholmod_sparse *cm_A21_sp = nullptr;
//        cholmod_dense *cm_A22_dn = nullptr;
//        cholmod_dense *cm_A12_dn = nullptr;
//        cholmod_factor *cm_L = nullptr;
//        cholmod_dense *cm_A11iA12_dn = nullptr;
//
//        setSymmetric(cm_A11_sp, A11_sp, false);
//        setSymmetric(cm_A21_sp, A21t_sp, false);
//        update(cm_A22_dn, A22t_dn);
//        update(cm_A12_dn, A12t_dn);
//
//        double alpha[2] = { -1, 0 };
//        double beta[2]  = {  1, 0 };
//
//        cm_L = _analyze<I>(cm_A11_sp, ext->cholmod.cm_common);
//        _factorize<I>(cm_L, cm_A11_sp, ext->cholmod.cm_common);
//        cm_A11iA12_dn = _solve<I>(CHOLMOD_A, cm_L, cm_A12_dn, ext->cholmod.cm_common);
//        _apply<I>(cm_A22_dn, cm_A21_sp, cm_A11iA12_dn, alpha, beta, ext->cholmod.cm_common);
//
//        if constexpr (std::is_same_v<T,double>) { sc.type = Matrix_Type::REAL_SYMMETRIC_POSITIVE_DEFINITE; }
//        if constexpr (std::is_same_v<T,std::complex<double>>) { sc.type = Matrix_Type::COMPLEX_HERMITIAN_POSITIVE_DEFINITE; }
//        sc.shape = Matrix_Shape::UPPER;
//        sc.resize(A22t_dn);
//        for(I r = 0, i = 0; r < sc.nrows; ++r) {
//            for(I c = r; c < sc.ncols; ++c, ++i) {
//                sc.vals[i] = A22t_dn.vals[r * sc.ncols + c];
//            }
//        }
//
//        delete cm_A11_sp;
//        delete cm_A21_sp;
//        delete cm_A22_dn;
//        delete cm_A12_dn;
//        _free<I>(cm_L, ext->cholmod.cm_common);
//        _free<I>(cm_A11iA12_dn, ext->cholmod.cm_common);
//    } break;
//
//    case Solver_External_Representation<T, I>::SOLVER::UMFPACK: {
//        eslog::error("implement DirectSparseSolver<T, I>::getSC\n");
//    } break;
//    }
}

}

#include "math/wrappers/math.spsolver.inst.hpp"

#endif
#endif
