
#ifndef SRC_WRAPPERS_SUITESPARSE_W_SUITESPARSE_CHOLMOD_H_
#define SRC_WRAPPERS_SUITESPARSE_W_SUITESPARSE_CHOLMOD_H_

#ifdef SUITESPARSE_HEADER_SUBDIR
#include "suitesparse/cholmod.h"
#else // SUITESPARSE_HEADER_DIRECT
#include "cholmod.h"
#endif

#include "math/math.h"

#include <cstring>

namespace espreso {

template <typename I> inline void _start(cholmod_common &common);
template <> inline void _start<int>(cholmod_common &common) { cholmod_start(&common); }
template <> inline void _start<long>(cholmod_common &common) { cholmod_l_start(&common); }

template <typename I> inline cholmod_sparse* _copy(cholmod_sparse* A, int stype, int mode, cholmod_common &common);
template <> inline cholmod_sparse* _copy<int>(cholmod_sparse* A, int stype, int mode, cholmod_common &common) { return cholmod_copy(A, stype, mode, &common); }
template <> inline cholmod_sparse* _copy<long>(cholmod_sparse* A, int stype, int mode, cholmod_common &common) { return cholmod_l_copy(A, stype, mode, &common); }

template <typename I> inline cholmod_factor* _analyze(cholmod_sparse* A, cholmod_common &common);
template <> inline cholmod_factor* _analyze<int>(cholmod_sparse* A, cholmod_common &common) { A->itype = CHOLMOD_INT; return cholmod_analyze(A, &common); }
template <> inline cholmod_factor* _analyze<long>(cholmod_sparse* A, cholmod_common &common) { A->itype = CHOLMOD_LONG; return cholmod_l_analyze(A, &common); }

template <typename I> inline void _factorize(cholmod_factor *L, cholmod_sparse *A, cholmod_common &common);
template <> inline void _factorize<int>(cholmod_factor *L, cholmod_sparse *A, cholmod_common &common) { cholmod_factorize(A, L, &common); }
template <> inline void _factorize<long>(cholmod_factor *L, cholmod_sparse *A, cholmod_common &common) { cholmod_l_factorize(A, L, &common); }

template <typename I> inline cholmod_dense* _solve(int sys, cholmod_factor *L, cholmod_dense *b, cholmod_common &common);
template <> inline cholmod_dense* _solve<int>(int sys, cholmod_factor *L, cholmod_dense *b, cholmod_common &common) { return cholmod_solve(sys, L, b, &common); }
template <> inline cholmod_dense* _solve<long>(int sys, cholmod_factor *L, cholmod_dense *b, cholmod_common &common) { return cholmod_l_solve(sys, L, b, &common); }

template <typename I> inline void _apply(cholmod_dense* Y, cholmod_sparse* A, cholmod_dense *X, double alpha[2], double beta[2], cholmod_common &common);
template <> inline void _apply<int>(cholmod_dense* Y, cholmod_sparse* A, cholmod_dense *X, double alpha[2], double beta[2], cholmod_common &common) { cholmod_sdmult(A, 0, alpha, beta, X, Y, &common); }
template <> inline void _apply<long>(cholmod_dense* Y, cholmod_sparse* A, cholmod_dense *X, double alpha[2], double beta[2], cholmod_common &common) { cholmod_l_sdmult(A, 0, alpha, beta, X, Y, &common); }

template <typename I> inline cholmod_sparse* _transpose(cholmod_sparse* A, cholmod_common &common);
template <> inline cholmod_sparse* _transpose<int>(cholmod_sparse* A, cholmod_common &common) { return cholmod_transpose(A, 1, &common); }
template <> inline cholmod_sparse* _transpose<long>(cholmod_sparse* A, cholmod_common &common) { return cholmod_l_transpose(A, 1, &common); }

template <typename I> inline cholmod_sparse* _multiply(cholmod_sparse* A, cholmod_sparse* B, cholmod_common &common);
template <> inline cholmod_sparse* _multiply<int>(cholmod_sparse* A, cholmod_sparse* B, cholmod_common &common) { return cholmod_ssmult(A, B, 1, true, true, &common); }
template <> inline cholmod_sparse* _multiply<long>(cholmod_sparse* A, cholmod_sparse* B, cholmod_common &common) { return cholmod_l_ssmult(A, B, 1, true, true, &common); }

template <typename I> inline cholmod_sparse* _factorToSparse(cholmod_factor *L, cholmod_common &common);
template <> inline cholmod_sparse* _factorToSparse<int>(cholmod_factor *L, cholmod_common &common) { return cholmod_factor_to_sparse(L, &common); }
template <> inline cholmod_sparse* _factorToSparse<long>(cholmod_factor *L, cholmod_common &common) { return cholmod_l_factor_to_sparse(L, &common); }

template <typename I> inline cholmod_factor* _copyFactor(cholmod_factor *in, cholmod_common &common);
template <> inline cholmod_factor* _copyFactor<int>(cholmod_factor *in, cholmod_common &common) { return cholmod_copy_factor(in, &common); }
template <> inline cholmod_factor* _copyFactor<long>(cholmod_factor *in, cholmod_common &common) { return cholmod_l_copy_factor(in, &common); }

template <typename I> inline void _changeFactor(int to_xtype, bool to_ll, bool to_super, bool to_packed, bool to_monotonic, cholmod_factor *L, cholmod_common &common);
template <> inline void _changeFactor<int>(int to_xtype, bool to_ll, bool to_super, bool to_packed, bool to_monotonic, cholmod_factor *L, cholmod_common &common) { cholmod_change_factor(to_xtype, to_ll, to_super, to_packed, to_monotonic, L, &common); }
template <> inline void _changeFactor<int64_t>(int to_xtype, bool to_ll, bool to_super, bool to_packed, bool to_monotonic, cholmod_factor *L, cholmod_common &common) { cholmod_l_change_factor(to_xtype, to_ll, to_super, to_packed, to_monotonic, L, &common); }

template <typename I> inline void _resymbol(cholmod_sparse *A, I *fset, size_t fsize, int pack, cholmod_factor *L, cholmod_common &common);
template <> inline void _resymbol<int>(cholmod_sparse *A, int *fset, size_t fsize, int pack, cholmod_factor *L, cholmod_common &common) { cholmod_resymbol(A, fset, fsize, pack, L, &common); }
template <> inline void _resymbol<int64_t>(cholmod_sparse *A, int64_t *fset, size_t fsize, int pack, cholmod_factor *L, cholmod_common &common) { cholmod_l_resymbol(A, fset, fsize, pack, L, &common); }

template <typename I> inline void _finish(cholmod_common &common);
template <> inline void _finish<int>(cholmod_common &common) { cholmod_finish(&common); }
template <> inline void _finish<long>(cholmod_common &common) { cholmod_l_finish(&common); }

template <typename I, typename O> inline void _free(O* &object, cholmod_common &common);
template <> inline void _free<int, cholmod_sparse>(cholmod_sparse* &object, cholmod_common &common) { cholmod_free_sparse(&object, &common); }
template <> inline void _free<long, cholmod_sparse>(cholmod_sparse* &object, cholmod_common &common) { cholmod_l_free_sparse(&object, &common); }
template <> inline void _free<int, cholmod_dense>(cholmod_dense* &object, cholmod_common &common) { cholmod_free_dense(&object, &common); }
template <> inline void _free<long, cholmod_dense>(cholmod_dense* &object, cholmod_common &common) { cholmod_l_free_dense(&object, &common); }
template <> inline void _free<int, cholmod_factor>(cholmod_factor* &object, cholmod_common &common) { cholmod_free_factor(&object, &common); }
template <> inline void _free<long, cholmod_factor>(cholmod_factor* &object, cholmod_common &common) { cholmod_l_free_factor(&object, &common); }

template <typename T> constexpr int _getCholmodXtype();
template <> constexpr int _getCholmodXtype<float>() { return CHOLMOD_REAL; }
template <> constexpr int _getCholmodXtype<double>() { return CHOLMOD_REAL; }
template <> constexpr int _getCholmodXtype<std::complex<float>>() { return CHOLMOD_COMPLEX; }
template <> constexpr int _getCholmodXtype<std::complex<double>>() { return CHOLMOD_COMPLEX; }

template <typename T> constexpr int _getCholmodDtype();
template <> constexpr int _getCholmodDtype<float>() { return CHOLMOD_SINGLE; }
template <> constexpr int _getCholmodDtype<double>() { return CHOLMOD_DOUBLE; }
template <> constexpr int _getCholmodDtype<std::complex<float>>() { return CHOLMOD_SINGLE; }
template <> constexpr int _getCholmodDtype<std::complex<double>>() { return CHOLMOD_DOUBLE; }

template<typename I> constexpr int _getCholmodItype();
template<> constexpr int _getCholmodItype<int32_t>() { return CHOLMOD_INT; }
template<> constexpr int _getCholmodItype<int64_t>() { return CHOLMOD_LONG; }

constexpr int _getCholmodStype(Matrix_Shape shape)
{
    switch(shape)
    {
    case Matrix_Shape::UPPER:
        return 1;
    case Matrix_Shape::LOWER:
        return -1;
    case Matrix_Shape::FULL:
        return 0;
    default:
        return 0;
    }
}

template <typename T, typename I>
static inline void setSymmetric(cholmod_sparse* &A, const Matrix_CSR<T, I> &M)
{
    if (!A) A = new cholmod_sparse();
    A->nrow = M.ncols; // CSR -> CSC, it works only for symmetric matrices
    A->ncol = M.nrows;
    A->nzmax = M.nnz;

    A->p = M.rows;
    A->i = M.cols;
    A->x = M.vals;

    A->stype = (-1) * _getCholmodStype(M.shape); // (-1)* <=> CSR -> CSC
    A->xtype = CHOLMOD_PATTERN;
    A->dtype = _getCholmodDtype<T>();
    A->sorted = 1;
    A->packed = 1;
}

template <typename T, typename I>
static inline void setAsymmetric(cholmod_sparse* &A, cholmod_common &common, const Matrix_CSR<T, I> &M)
{
    if (A) delete A;
    cholmod_sparse* At = new cholmod_sparse();
    At->nrow = M.ncols;
    At->ncol = M.nrows;
    At->nzmax = M.nnz;

    At->p = M.rows;
    At->i = M.cols;
    At->x = M.vals;

    At->stype = (-1) * _getCholmodStype(M.shape); // (-1)* <=> CSR -> CSC
    At->xtype = _getCholmodXtype<T>();
    At->dtype = _getCholmodDtype<T>();
    At->sorted = 1;
    At->packed = 1;

    // CSR -> CSC
    At = _transpose<esint>(A, common);
    delete At;
}

template <typename T, typename I>
static inline void updateSymmetric(cholmod_sparse *A, const Matrix_CSR<T, I> &M)
{
    A->x = M.vals;
    A->xtype = _getCholmodXtype<T>();
}

template <typename T, typename I>
static inline void update(cholmod_dense* &A, const Matrix_Dense<T, I> &M)
{
    if (!A) A = new cholmod_dense();
    A->nrow = M.ncols; // row-major -> column-major, therefore transposed. JUST array-transposed. NOT conjugate-transposed
    A->ncol = M.nrows;
    A->nzmax = M.ncols * M.nrows;
    A->d = M.ncols;
    A->x = M.vals;
    A->xtype = _getCholmodXtype<T>();
    A->dtype = _getCholmodDtype<T>();
}

template <typename T, typename I>
static inline void update(cholmod_dense* &A, const Vector_Dense<T, I> &v)
{
    if (!A) A = new cholmod_dense();
    A->nrow = v.size;
    A->ncol = 1;
    A->nzmax = v.size;
    A->d = v.size;
    A->x = v.vals;
    A->xtype = _getCholmodXtype<T>();
    A->dtype = _getCholmodDtype<T>();
}

template <typename T, typename I>
static inline void _extractUpper(cholmod_sparse* &A, cholmod_common &common, Matrix_CSR<T, I> &M)
{
    cholmod_sparse* upA = _copy<esint>(A, -1, 1, common); // CSC -> CSR
    M.shape = Matrix_Shape::UPPER;
    M.type = Matrix_Type::REAL_SYMMETRIC_INDEFINITE;
    M.resize(upA->ncol, upA->nrow, upA->nzmax);
    std::copy((esint*)upA->p, (esint*)upA->p + M.nrows + 1, M.rows);
    std::copy((esint*)upA->i, (esint*)upA->i + M.nnz, M.cols);
    std::copy((T*)upA->x, (T*)upA->x + M.nnz, M.vals);
    delete upA;
}

template <typename T, typename I>
static inline void extract(cholmod_dense *A, cholmod_common &common, Matrix_Dense<T, I> &M)
{
    memcpy(M.vals, A->x, sizeof(T) * M.nrows * M.ncols);
    _free<esint>(A, common);
}

template <typename T, typename I>
static inline void extract(cholmod_dense *A, cholmod_common &common, Vector_Dense<T, I> &v)
{
    memcpy(v.vals, A->x, sizeof(T) * v.size);
    _free<esint>(A, common);
}

}

#endif /* SRC_WRAPPERS_SUITESPARSE_W_SUITESPARSE_CHOLMOD_H_ */
