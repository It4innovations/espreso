
#include "math/math.h"
#include "esinfo/eslog.h"
<<<<<<< HEAD
#include "gpu/gpu_management.h"
=======
>>>>>>> 8a13c706 (ENH: new spsolver interface, new suitesparse spsolver implementation, other changes to make it work)
#include "basis/utilities/utils.h"

#include <complex>

#ifndef HAVE_MKL
#ifdef HAVE_SUITESPARSE

<<<<<<< HEAD
#include "math/wrappers/math.spsolver.h"
=======
>>>>>>> 8a13c706 (ENH: new spsolver interface, new suitesparse spsolver implementation, other changes to make it work)
#include "w.suitesparse.cholmod.h"

namespace espreso {

template<typename T, typename I>
struct Solver_External_Representation {
<<<<<<< HEAD
	cholmod_common cm_common;
	cholmod_factor * cm_factor_super = nullptr;
	cholmod_factor * cm_factor_simpl = nullptr;
	cholmod_sparse * cm_matrix_view = nullptr;
    const Matrix_CSR<T, I> * matrix = nullptr;
	Vector_Dense<I> map_simpl_super;
	char zerodrop;
	int stage = 0; // 0 = completely uninitialized, 1 = initialized without matrix, 2 = have matrix, 3 = symbolic factorization done, 4 = numeric factorization done
};


=======
	char zerodrop;
	int stage = 0; // 0 = completely uninitialized, 1 = initialized without matrix, 2 = matrix set but not factorized, 3 = symbolic factorization done, 4 = numeric factorization done
	cholmod_common * cm_common = nullptr;
	cholmod_factor * cm_factor_super = nullptr;
	cholmod_factor * cm_factor_simpl = nullptr;
	cholmod_sparse * cm_matrix_view = nullptr;
	Vector_Dense<I> map_simpl_super;
};




template<typename I>
static void my_cholmod_start(cholmod_common * common)
{
    if constexpr(std::is_same_v<I,int32_t>) cholmod_start(common);
    if constexpr(std::is_same_v<I,int64_t>) cholmod_l_start(common);
}
template<typename I>
static cholmod_factor * my_cholmod_analyze(cholmod_sparse * A, cholmod_common * common)
{
    if constexpr(std::is_same_v<I,int32_t>) return cholmod_analyze(A, common);
    if constexpr(std::is_same_v<I,int64_t>) return cholmod_l_analyze(A, common);
}
template<typename I>
static void my_cholmod_factorize(cholmod_sparse * A, cholmod_factor * F, cholmod_common * common)
{
    if constexpr(std::is_same_v<I,int32_t>) cholmod_factorize(A, F, common);
    if constexpr(std::is_same_v<I,int64_t>) cholmod_l_factorize(A, F, common);
}
template<typename I>
static cholmod_sparse * my_cholmod_factor_to_sparse(cholmod_factor * F, cholmod_common * common)
{
    if constexpr(std::is_same_v<I,int32_t>) return cholmod_factor_to_sparse(F, common);
    if constexpr(std::is_same_v<I,int64_t>) return cholmod_l_factor_to_sparse(F, common);
}
template<typename I>
static cholmod_dense * my_cholmod_solve(int sys, cholmod_factor * L, cholmod_dense * B, cholmod_common * common)
{
    if constexpr(std::is_same_v<I,int32_t>) return cholmod_solve(sys, L, B, common);
    if constexpr(std::is_same_v<I,int64_t>) return cholmod_l_solve(sys, L, B, common);
}
template<typename I>
static void my_cholmod_free_sparse(cholmod_sparse ** A, cholmod_common * common)
{
    if constexpr(std::is_same_v<I,int32_t>) cholmod_free_sparse(A, common);
    if constexpr(std::is_same_v<I,int64_t>) cholmod_l_free_sparse(A, common);
}
template<typename I>
static void my_cholmod_free_factor(cholmod_factor ** F, cholmod_common * common)
{
    if constexpr(std::is_same_v<I,int32_t>) cholmod_free_factor(F, common);
    if constexpr(std::is_same_v<I,int64_t>) cholmod_l_free_factor(F, common);
}
template<typename I>
static int my_cholmod_sdmult(cholmod_sparse * A, int transpose, double alpha[2], double beta[2], cholmod_dense * X, cholmod_dense * Y, cholmod_common * Common)
{
    if constexpr(std::is_same_v<I,int32_t>) return cholmod_sdmult(A, transpose, alpha, beta, X, Y, Common);
    if constexpr(std::is_same_v<I,int64_t>) return cholmod_l_sdmult(A, transpose, alpha, beta, X, Y, Common);
}
template<typename I>
static cholmod_dense * my_cholmod_sparse_to_dense(cholmod_sparse * A, cholmod_common * Common)
{
    if constexpr(std::is_same_v<I,int32_t>) return cholmod_sparse_to_dense(A, Common);
    if constexpr(std::is_same_v<I,int64_t>) return cholmod_l_sparse_to_dense(A, Common);
}
template<typename I>
static cholmod_dense * my_cholmod_allocate_dense(size_t nrow, size_t ncol, size_t d, int xtype, cholmod_common *Common)
{
    if constexpr(std::is_same_v<I,int32_t>) return cholmod_allocate_dense(nrow, ncol, d, xtype, Common);
    if constexpr(std::is_same_v<I,int64_t>) return cholmod_l_allocate_dense(nrow, ncol, d, xtype, Common);
}
template<typename I>
static int my_cholmod_free_dense(cholmod_dense **X, cholmod_common *Common)
{
    if constexpr(std::is_same_v<I,int32_t>) return cholmod_free_dense(X, Common);
    if constexpr(std::is_same_v<I,int64_t>) return cholmod_l_free_dense(X, Common);
}
template<typename I>
static cholmod_factor * my_cholmod_copy_factor(cholmod_factor * F, cholmod_common * Common)
{
    if constexpr(std::is_same_v<I,int32_t>) return cholmod_copy_factor(F, Common);
    if constexpr(std::is_same_v<I,int64_t>) return cholmod_l_copy_factor(F, Common);
}
template<typename I>
static int my_cholmod_change_factor(int to_xtype, int to_ll, int to_super, int to_packed, int to_monotonic, cholmod_factor *L, cholmod_common *Common)
{
    if constexpr(std::is_same_v<I,int32_t>) return cholmod_change_factor(to_xtype, to_ll, to_super, to_packed, to_monotonic, L, Common);
    if constexpr(std::is_same_v<I,int64_t>) return cholmod_l_change_factor(to_xtype, to_ll, to_super, to_packed, to_monotonic, L, Common);
}
template<typename I>
static int my_cholmod_resymbol(cholmod_sparse *A, I *fset, size_t fsize, int pack, cholmod_factor *L, cholmod_common *Common)
{
    if constexpr(std::is_same_v<I,int32_t>) return cholmod_resymbol(A, fset, fsize, pack, L, Common);
    if constexpr(std::is_same_v<I,int64_t>) return cholmod_l_resymbol(A, fset, fsize, pack, L, Common);
}
template<typename T>
static constexpr int my_cholmod_dtype()
{
    if constexpr(std::is_same_v<T,double>)               return CHOLMOD_DOUBLE;
    if constexpr(std::is_same_v<T,std::complex<double>>) return CHOLMOD_DOUBLE;
    if constexpr(std::is_same_v<T,float>)                return CHOLMOD_SINGLE;
    if constexpr(std::is_same_v<T,std::complex<float>>)  return CHOLMOD_SINGLE;
}
template<typename T>
static constexpr int my_cholmod_xtype()
{
    if constexpr(std::is_same_v<T,std::complex<float>>)  return CHOLMOD_COMPLEX;
    if constexpr(std::is_same_v<T,std::complex<double>>) return CHOLMOD_COMPLEX;
    if constexpr(std::is_same_v<T,float>)                return CHOLMOD_REAL;
    if constexpr(std::is_same_v<T,double>)               return CHOLMOD_REAL;
}
template<typename I>
static constexpr int my_cholmod_itype()
{
    if constexpr(std::is_same_v<I,int32_t>) return CHOLMOD_INT;
    if constexpr(std::is_same_v<I,int64_t>) return CHOLMOD_LONG;
}



>>>>>>> 8a13c706 (ENH: new spsolver interface, new suitesparse spsolver implementation, other changes to make it work)
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
	return Solver_Factors::SYMMETRIC_UPPER;
}

template <typename T, typename I>
DirectSparseSolver<T, I>::DirectSparseSolver()
{
	ext = std::make_unique<Solver_External_Representation<T,I>>();

    ext->zerodrop = 'D'; // Drop, Keep

<<<<<<< HEAD
    _start<I>(ext->cm_common);
    ext->cm_common.final_ll = 1;
    ext->cm_common.nthreads_max = 1;
    ext->cm_common.nmethods = 1;
    ext->cm_common.method[0].ordering = CHOLMOD_METIS;
    ext->cm_common.itype = _getCholmodItype<I>();
    ext->cm_common.supernodal = CHOLMOD_SUPERNODAL;
=======
    ext->cm_common = new cholmod_common();
    my_cholmod_start<I>(ext->cm_common);
    ext->cm_common->final_ll = 1;
    ext->cm_common->nthreads_max = 1;
    ext->cm_common->nmethods = 1;
    ext->cm_common->method[0].ordering = CHOLMOD_METIS;
    ext->cm_common->itype = my_cholmod_itype<I>();
    ext->cm_common->supernodal = CHOLMOD_SUPERNODAL;
>>>>>>> 8a13c706 (ENH: new spsolver interface, new suitesparse spsolver implementation, other changes to make it work)

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

<<<<<<< HEAD
    if(ext->cm_factor_super != nullptr) _free<I>(ext->cm_factor_super, ext->cm_common);

    if(ext->cm_factor_super != nullptr) _free<I>(ext->cm_factor_simpl, ext->cm_common);

    _finish<I>(ext->cm_common);
=======
    if(ext->cm_factor_super != nullptr) my_cholmod_free_factor<I>(&ext->cm_factor_super, ext->cm_common);

    if(ext->cm_factor_super != nullptr) my_cholmod_free_factor<I>(&ext->cm_factor_simpl, ext->cm_common);

    if(ext->cm_common != nullptr) cholmod_finish(ext->cm_common);
    delete ext->cm_common;
>>>>>>> 8a13c706 (ENH: new spsolver interface, new suitesparse spsolver implementation, other changes to make it work)
}

template <typename T, typename I>
void DirectSparseSolver<T, I>::commit(const Matrix_CSR<T,I> &a)
{
	if(a.nrows != a.ncols) eslog::error("commit: matrix has to be square\n");
	if(a.type != Matrix_Type::REAL_SYMMETRIC_POSITIVE_DEFINITE && a.type != Matrix_Type::COMPLEX_HERMITIAN_POSITIVE_DEFINITE) eslog::error("commit: matrix has to be SPD or HPD\n");
	if(a.shape != Matrix_Shape::UPPER) eslog::error("commit: CSR matrix has to be upper triangular\n");

<<<<<<< HEAD
    if(ext->stage < 1) eslog::error("commit: invalid order of operations in spsolver\n");

    ext->matrix = &a;
=======
    if(ext->stage < 1) eslog::error("commit: invalid order of operations in solver\n");
>>>>>>> 8a13c706 (ENH: new spsolver interface, new suitesparse spsolver implementation, other changes to make it work)

	if(ext->cm_matrix_view == nullptr)
    {
        ext->cm_matrix_view = new cholmod_sparse();
        ext->cm_matrix_view->nrow = a.ncols;
        ext->cm_matrix_view->ncol = a.nrows;
        ext->cm_matrix_view->nzmax = a.nnz;
        ext->cm_matrix_view->nz = nullptr;
        ext->cm_matrix_view->z = nullptr;
<<<<<<< HEAD
        ext->cm_matrix_view->stype = -1; // UPPER in CSR, but LOWER in CSC
        ext->cm_matrix_view->itype = _getCholmodItype<I>();
        ext->cm_matrix_view->xtype = _getCholmodXtype<T>();
        ext->cm_matrix_view->dtype = _getCholmodDtype<T>();
=======
        ext->cm_matrix_view->stype = -1;
        ext->cm_matrix_view->itype = my_cholmod_itype<I>();
        ext->cm_matrix_view->xtype = my_cholmod_xtype<T>();
        ext->cm_matrix_view->dtype = my_cholmod_dtype<T>();
>>>>>>> 8a13c706 (ENH: new spsolver interface, new suitesparse spsolver implementation, other changes to make it work)
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

<<<<<<< HEAD
    if(ext->stage != 2) throw std::runtime_error("symbolicFactorization: invalid order of operations in spsolver\n");
	
	ext->cm_factor_super = _analyze<I>(ext->cm_matrix_view, ext->cm_common);

    if(ext->cm_factor_super->xsize > utils::get_max_val_no_precision_loss_in_fp<T>()) eslog::error("symbolicFactorization: factor nnz too large for my super->simpl map\n");
    ext->cm_factor_simpl = _copyFactor<I>(ext->cm_factor_super, ext->cm_common);
    _changeFactor<I>(_getCholmodXtype<T>(), true, true, true, true, ext->cm_factor_simpl, ext->cm_common);
    for(size_t i = 0; i < ext->cm_factor_simpl->xsize; i++) reinterpret_cast<T*>(ext->cm_factor_simpl->x)[i] = static_cast<T>(i);
    _changeFactor<I>(_getCholmodXtype<T>(), true, false, true, true, ext->cm_factor_simpl, ext->cm_common);
    if(ext->zerodrop == 'D') _resymbol<I>(ext->cm_matrix_view, nullptr, 0, 1, ext->cm_factor_simpl, ext->cm_common);
=======
    if(ext->stage != 2) throw std::runtime_error("symbolicFactorization: invalid order of operations in solver\n");
	
	ext->cm_factor_super = my_cholmod_analyze<I>(ext->cm_matrix_view, ext->cm_common);

    if(ext->cm_factor_super->xsize > utils::get_max_val_no_precision_loss_in_fp<T>()) eslog::error("symbolicFactorization: factor nnz too large for my super->simpl map\n");
    ext->cm_factor_simpl = my_cholmod_copy_factor<I>(ext->cm_factor_super, ext->cm_common);
    my_cholmod_change_factor<I>(CHOLMOD_REAL, 1, 1, 1, 1, ext->cm_factor_simpl, ext->cm_common);
    for(size_t i = 0; i < ext->cm_factor_simpl->xsize; i++) reinterpret_cast<T*>(ext->cm_factor_simpl->x)[i] = static_cast<T>(i);
    my_cholmod_change_factor<I>(CHOLMOD_REAL, 1, 0, 1, 1, ext->cm_factor_simpl, ext->cm_common);
    if(ext->zerodrop == 'D') my_cholmod_resymbol<I>(ext->cm_matrix_view, nullptr, 0, 1, ext->cm_factor_simpl, ext->cm_common);
>>>>>>> 8a13c706 (ENH: new spsolver interface, new suitesparse spsolver implementation, other changes to make it work)
    ext->map_simpl_super.resize(ext->cm_factor_simpl->nzmax);
    for(I i = 0; i < ext->map_simpl_super.size; i++) ext->map_simpl_super.vals[i] = static_cast<I>(std::real(reinterpret_cast<T*>(ext->cm_factor_simpl->x)[i]));
    
    ext->stage = 3;
}

template <typename T, typename I>
void DirectSparseSolver<T, I>::numericalFactorization()
{
<<<<<<< HEAD
	if(ext->stage < 3) eslog::error("numericalFactorization: invalid order of operations in spsolver\n");

    _factorize<I>(ext->cm_factor_super, ext->cm_matrix_view, ext->cm_common);
=======
	if(ext->stage < 3) eslog::error("numericalFactorization: invalid order of operations in solver\n");

    my_cholmod_factorize<I>(ext->cm_matrix_view, ext->cm_factor_super, ext->cm_common);
>>>>>>> 8a13c706 (ENH: new spsolver interface, new suitesparse spsolver implementation, other changes to make it work)

    ext->stage = 4;
}

template <typename T, typename I>
void DirectSparseSolver<T, I>::solve(Vector_Dense<T, I> &rhs, Vector_Dense<T, I> &solution, int sparsity)
{
<<<<<<< HEAD
	if(ext->stage < 4) eslog::error("solve: invalid order of operations in spsolver\n");

=======
>>>>>>> 8a13c706 (ENH: new spsolver interface, new suitesparse spsolver implementation, other changes to make it work)
	cholmod_dense cm_rhs;
	cm_rhs.nrow = rhs.size;
	cm_rhs.ncol = 1;
	cm_rhs.d = rhs.size;
	cm_rhs.nzmax = rhs.size;
	cm_rhs.x = rhs.vals;
<<<<<<< HEAD
	cm_rhs.xtype = _getCholmodXtype<T>();
	cm_rhs.dtype = _getCholmodDtype<T>();

	cholmod_dense * cm_sol = _solve<I>(CHOLMOD_A, ext->cm_factor_super, &cm_rhs, ext->cm_common);
=======
	cm_rhs.xtype = my_cholmod_xtype<T>();
	cm_rhs.dtype = my_cholmod_dtype<T>();

	cholmod_dense * cm_sol = my_cholmod_solve<I>(CHOLMOD_A, ext->cm_factor_super, &cm_rhs, ext->cm_common);
>>>>>>> 8a13c706 (ENH: new spsolver interface, new suitesparse spsolver implementation, other changes to make it work)

	solution.resize(cm_sol->nrow);
	std::copy_n(reinterpret_cast<T*>(cm_sol->x), cm_sol->nrow, solution.vals);

<<<<<<< HEAD
	_free<I>(cm_sol, ext->cm_common);
=======
	my_cholmod_free_dense<I>(&cm_sol, ext->cm_common);
>>>>>>> 8a13c706 (ENH: new spsolver interface, new suitesparse spsolver implementation, other changes to make it work)
}

template <typename T, typename I>
void DirectSparseSolver<T, I>::solve(Matrix_Dense<T, I> &rhs, Matrix_Dense<T, I> &solution, int sparsity)
{
<<<<<<< HEAD
	if(ext->stage < 4) eslog::error("solve: invalid order of operations in spsolver\n");

=======
>>>>>>> 8a13c706 (ENH: new spsolver interface, new suitesparse spsolver implementation, other changes to make it work)
	cholmod_dense cm_rhs;
	cm_rhs.nrow = rhs.ncols;
	cm_rhs.ncol = rhs.nrows;
	cm_rhs.d = rhs.get_ld();
	cm_rhs.nzmax = cm_rhs.d * rhs.nrows;
	cm_rhs.x = rhs.vals;
<<<<<<< HEAD
	cm_rhs.xtype = _getCholmodXtype<T>();
	cm_rhs.dtype = _getCholmodDtype<T>();

	cholmod_dense * cm_sol = _solve<I>(CHOLMOD_A, ext->cm_factor_super, &cm_rhs, ext->cm_common);
=======
	cm_rhs.xtype = my_cholmod_xtype<T>();
	cm_rhs.dtype = my_cholmod_dtype<T>();

	cholmod_dense * cm_sol = my_cholmod_solve<I>(CHOLMOD_A, ext->cm_factor_super, &cm_rhs, ext->cm_common);
>>>>>>> 8a13c706 (ENH: new spsolver interface, new suitesparse spsolver implementation, other changes to make it work)

	solution.resize(cm_sol->ncol, cm_sol->d);
	solution.ncols = cm_sol->nrow;
	std::copy_n(reinterpret_cast<T*>(cm_sol->x), cm_sol->nzmax, solution.vals);

<<<<<<< HEAD
	_free<I>(cm_sol, ext->cm_common);
=======
	my_cholmod_free_dense<I>(&cm_sol, ext->cm_common);
>>>>>>> 8a13c706 (ENH: new spsolver interface, new suitesparse spsolver implementation, other changes to make it work)
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
<<<<<<< HEAD
	if(ext->stage < 3) eslog::error("getFactorNnz: invalid order of operations in spsolver\n");
=======
	if(ext->stage < 3) eslog::error("getFactorNnz: invalid order of operations in solver\n");
>>>>>>> 8a13c706 (ENH: new spsolver interface, new suitesparse spsolver implementation, other changes to make it work)

    // https://github.com/DrTimothyAldenDavis/SuiteSparse/issues/523

    return ext->cm_factor_simpl->nzmax;
}

template <typename T, typename I>
<<<<<<< HEAD
template <typename A>
inline void DirectSparseSolver<T, I>::getFactorL(Matrix_CSR<T,I,A> &/*L*/, bool /*copyPattern*/, bool /*copyValues*/)
=======
void DirectSparseSolver<T, I>::getFactorL(Matrix_CSR<T,I> &/*L*/, bool /*copyPattern*/, bool /*copyValues*/)
>>>>>>> 8a13c706 (ENH: new spsolver interface, new suitesparse spsolver implementation, other changes to make it work)
{
	eslog::error("L factor is not provided\n");
}

template <typename T, typename I>
<<<<<<< HEAD
template <typename A>
inline void DirectSparseSolver<T, I>::getFactorU(Matrix_CSR<T,I,A> &U, bool copyPattern, bool copyValues)
{
	if(ext->stage < 3) eslog::error("getFactorU: invalid order of operations in spsolver\n");
    if(copyValues && ext->stage < 4) eslog::error("getFactorU: invalid order of operations in spsolver\n");
=======
void DirectSparseSolver<T, I>::getFactorU(Matrix_CSR<T,I> &U, bool copyPattern, bool copyValues)
{
	if(ext->stage < 3) eslog::error("getFactorU: invalid order of operations in solver\n");
    if(copyValues && ext->stage < 4) eslog::error("getFactorU: invalid order of operations in solver\n");
>>>>>>> 8a13c706 (ENH: new spsolver interface, new suitesparse spsolver implementation, other changes to make it work)
    if((size_t)U.nrows != ext->cm_factor_simpl->n || (size_t)U.ncols != ext->cm_factor_simpl->n) eslog::error("getFactorU: output matrix has wrong dimensions\n");

	U.resize(ext->cm_factor_simpl->n, ext->cm_factor_simpl->n, ext->cm_factor_simpl->nzmax);

    if(copyPattern) std::copy_n(static_cast<I*>(ext->cm_factor_simpl->p), ext->cm_factor_simpl->n+1, U.rows);
    if(copyPattern) std::copy_n(static_cast<I*>(ext->cm_factor_simpl->i), ext->cm_factor_simpl->nzmax, U.cols);
    if(copyValues) for(I i = 0; i < ext->map_simpl_super.size; i++) U.vals[i] = reinterpret_cast<T*>(ext->cm_factor_super->x)[ext->map_simpl_super.vals[i]];
}

template <typename T, typename I>
void DirectSparseSolver<T, I>::getPermutation(Permutation<I> &perm)
{
<<<<<<< HEAD
	if(ext->stage < 3) eslog::error("getPermutation: invalid order of operations in spsolver\n");
=======
	if(ext->stage < 3) eslog::error("getPermutation: invalid order of operations in solver\n");
>>>>>>> 8a13c706 (ENH: new spsolver interface, new suitesparse spsolver implementation, other changes to make it work)

	perm.resize(ext->cm_factor_simpl->n);

    std::copy_n(static_cast<I*>(ext->cm_factor_simpl->Perm), ext->cm_factor_simpl->n, perm.dst_to_src);

    if(ext->cm_factor_simpl->IPerm != nullptr) std::copy_n(static_cast<I*>(ext->cm_factor_simpl->IPerm), ext->cm_factor_simpl->n, perm.dst_to_src);
    else perm.invert(perm.dst_to_src, perm.src_to_dst);
}

template <typename T, typename I>
<<<<<<< HEAD
void DirectSparseSolver<T, I>::getPermutation(Vector_Dense<I> &perm)
{
	if(ext->stage < 3) eslog::error("getPermutation: invalid order of operations in spsolver\n");

	perm.resize(ext->cm_factor_simpl->n);

    std::copy_n(static_cast<I*>(ext->cm_factor_simpl->Perm), ext->cm_factor_simpl->n, perm.vals);
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

    cm_L = _analyze<esint>(cm_A11_sp, ext->cm_common);
    _factorize<esint>(cm_L, cm_A11_sp, ext->cm_common);
    cm_A11iA12_dn = _solve<esint>(CHOLMOD_A, cm_L, cm_A12_dn, ext->cm_common);
    _apply<esint>(cm_A22_dn, cm_A21_sp, cm_A11iA12_dn, alpha, beta, ext->cm_common);

    if constexpr (std::is_same_v<T,double>) { sc.type = Matrix_Type::REAL_SYMMETRIC_POSITIVE_DEFINITE; }
    if constexpr (std::is_same_v<T,std::complex<double>>) { sc.type = Matrix_Type::COMPLEX_HERMITIAN_POSITIVE_DEFINITE; }
    sc.shape = Matrix_Shape::UPPER;
    sc.resize(A22t_dn);
    for(esint r = 0, i = 0; r < sc.nrows; ++r) {
        for(esint c = r; c < sc.ncols; ++c, ++i) {
            sc.vals[i] = A22t_dn.vals[r * sc.ncols + c];
        }
    }

    delete cm_A11_sp;
    delete cm_A21_sp;
    delete cm_A22_dn;
    delete cm_A12_dn;
    _free<esint>(cm_L, ext->cm_common);
    _free<esint>(cm_A11iA12_dn, ext->cm_common);
}

template struct DirectSparseSolver<float, int32_t>;
template struct DirectSparseSolver<double, int32_t>;
template struct DirectSparseSolver<float, int64_t>;
template struct DirectSparseSolver<double, int64_t>;
template struct DirectSparseSolver<std::complex<float>, int32_t>;
template struct DirectSparseSolver<std::complex<double>, int32_t>;
template struct DirectSparseSolver<std::complex<float>, int64_t>;
template struct DirectSparseSolver<std::complex<double>, int64_t>;

template void DirectSparseSolver<double, int32_t>::getFactorL<gpu::mgm::Ah>(Matrix_CSR<double, int32_t, gpu::mgm::Ah> &, bool, bool);
template void DirectSparseSolver<double, int32_t>::getFactorU<gpu::mgm::Ah>(Matrix_CSR<double, int32_t, gpu::mgm::Ah> &, bool, bool);
=======
void DirectSparseSolver<T, I>::getSC(Matrix_Dense<T,I> &sc)
{
    eslog::error("getSC: not implemented, todo\n");
    // old implementation, todo
	// computes the schur complement S = A22 - A21 * A11^{-1} * A12, where A = [A11, A12; A21, A22]

    // esint size_sc = sc.nrows;
    // esint size = matrix->nrows;
    // esint size_A11 = size - size_sc;

    // Matrix_CSR<T, I> A11_sp;
    // Matrix_CSR<T, I> A21t_sp; // = A12c_sp
    // Matrix_Dense<T, I> A22t_dn;
    // Matrix_Dense<T, I> A12t_dn;
    // SpBLAS<Matrix_CSR, T, I>::submatrix(*matrix, A11_sp, 0, size_A11, 0, size_A11);
    // SpBLAS<Matrix_CSR, T, I>::submatrix(*matrix, A21t_sp, 0, size_A11, size_A11, size, false, true); // = A12c_sp
    // SpBLAS<Matrix_CSR, T, I>::submatrix(*matrix, A22t_dn, size_A11, size, size_A11, size, true, false, true);
    // SpBLAS<Matrix_CSR, T, I>::submatrix(*matrix, A12t_dn, 0, size_A11, size_A11, size, true, false, true);

    // cholmod_common &cm_common = _solver->cholmod.common;
    // cholmod_sparse *cm_A11_sp = new cholmod_sparse();
    // cholmod_sparse *cm_A21_sp = new cholmod_sparse();
    // cholmod_dense *cm_A22_dn = new cholmod_dense();
    // cholmod_dense *cm_A12_dn = new cholmod_dense();
    // cholmod_factor *cm_L;
    // cholmod_dense *cm_A11iA12_dn;

    // setSymmetric(cm_A11_sp, A11_sp);
    // updateSymmetric(cm_A11_sp, A11_sp);
    // setSymmetric(cm_A21_sp, A21t_sp);
    // updateSymmetric(cm_A21_sp, A21t_sp);
    // update(cm_A22_dn, A22t_dn);
    // update(cm_A12_dn, A12t_dn);

    // double alpha[2] = {-1,0};
    // double beta[2] = {1,0};

    // _analyze<esint>(cm_L, cm_A11_sp, cm_common);
    // _factorize<esint>(cm_L, cm_A11_sp, cm_common);
    // _solve<esint>(cm_A11iA12_dn, cm_L, cm_A12_dn, cm_common);
    // _apply<esint>(cm_A22_dn, cm_A21_sp, cm_A11iA12_dn, alpha, beta, cm_common);

    // if constexpr (std::is_same_v<T,double>) { sc.type = Matrix_Type::REAL_SYMMETRIC_POSITIVE_DEFINITE; }
    // if constexpr (std::is_same_v<T,std::complex<double>>) { sc.type = Matrix_Type::COMPLEX_HERMITIAN_POSITIVE_DEFINITE; }
    // sc.shape = Matrix_Shape::UPPER;
    // sc.resize(A22t_dn);
    // for(esint r = 0, i = 0; r < sc.nrows; ++r) {
    //     for(esint c = r; c < sc.ncols; ++c, ++i) {
    //         sc.vals[i] = A22t_dn.vals[r * sc.ncols + c];
    //     }
    // }

    // delete cm_A11_sp;
    // delete cm_A21_sp;
    // delete cm_A22_dn;
    // delete cm_A12_dn;
    // _free<esint>(cm_L, cm_common);
    // _free<esint>(cm_A11iA12_dn, cm_common);
    // break;
}

template struct DirectSparseSolver<double, int>;
template struct DirectSparseSolver<std::complex<double>, int>;
>>>>>>> 8a13c706 (ENH: new spsolver interface, new suitesparse spsolver implementation, other changes to make it work)

}

#endif
#endif

