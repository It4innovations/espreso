
#include "math/wrappers/math.spblas.h"
#include "esinfo/eslog.h"
#include "math/operations/copy_dense.h"

#include <complex>

#ifdef HAVE_SUITESPARSE
#ifdef ESPRESO_USE_WRAPPER_SPBLAS_SUITESPARSE

#include "w.suitesparse.cholmod.h"

namespace espreso {

struct Matrix_SpBLAS_External_Representation {
    cholmod_sparse *A;
    cholmod_dense *X, *Y;
    double alpha[2] = { 1., 1. }, beta[2] = { 0., 0. };
    cholmod_common common;

    Matrix_SpBLAS_External_Representation(): A(new cholmod_sparse()), X(new cholmod_dense()), Y(new cholmod_dense()) { }
    ~Matrix_SpBLAS_External_Representation() { delete A; delete X; delete Y;}
};

template <template <typename, typename, typename> class Matrix, typename T, typename I>
SpBLAS<Matrix, T, I>::SpBLAS()
: matrix{}, _spblas{}
{

}

template <template <typename, typename, typename> class Matrix, typename T, typename I>
SpBLAS<Matrix, T, I>::~SpBLAS()
{
    if (_spblas) {
        _finish<I>(_spblas->common);
        delete _spblas;
    }
}

template <template <typename, typename, typename> class Matrix, typename T, typename I>
SpBLAS<Matrix, T, I>::SpBLAS(MatrixType &a, bool trans)
: matrix(&a), _spblas(nullptr)
{
    insert(a, trans);
}

template <template <typename, typename, typename> class Matrix, typename T, typename I>
void SpBLAS<Matrix, T, I>::insert(MatrixType &a, bool trans)
{
    matrix = &a;
    if (_spblas) {
        _finish<I>(_spblas->common);
        delete _spblas;
    }

    _spblas = new Matrix_SpBLAS_External_Representation();
    _start<I>(_spblas->common);
    if (a.shape == Matrix_Shape::FULL) {
        setAsymmetric(_spblas->A, _spblas->common, *matrix, trans);
        matrix = nullptr;
    } else {
        setSymmetric(_spblas->A, *matrix, trans);
    }
}

template <>
void SpBLAS<Matrix_CSR, float, int>::apply(Vector_Dense<float> &y, const float &alpha, const float &beta, const Vector_Dense<float> &x)
{
    if (_spblas->A->nrow == 0 || _spblas->A->ncol == 0) return;

    update(_spblas->X, x);
    update(_spblas->Y, y);
    _spblas->alpha[0] = alpha;
    _spblas->beta[0] = beta;
    _apply<int>(_spblas->Y, _spblas->A, _spblas->X, _spblas->alpha, _spblas->beta, _spblas->common);
}

template <>
void SpBLAS<Matrix_CSR, double, int>::apply(Vector_Dense<double> &y, const double &alpha, const double &beta, const Vector_Dense<double> &x)
{
    if (_spblas->A->nrow == 0 || _spblas->A->ncol == 0) return;

    update(_spblas->X, x);
    update(_spblas->Y, y);
    _spblas->alpha[0] = alpha;
    _spblas->beta[0] = beta;
    _apply<int>(_spblas->Y, _spblas->A, _spblas->X, _spblas->alpha, _spblas->beta, _spblas->common);
}

template <>
void SpBLAS<Matrix_CSR, std::complex<double>, int>::apply(Vector_Dense<std::complex<double>, int> &y, const std::complex<double> &alpha, const std::complex<double> &beta, const Vector_Dense<std::complex<double> > &x)
{
    if (_spblas->A->nrow == 0 || _spblas->A->ncol == 0) return;

    update(_spblas->X, x);
    update(_spblas->Y, y);
    _spblas->alpha[0] = alpha.real();
    _spblas->alpha[1] = alpha.imag();
    _spblas->beta[0] = beta.real();
    _spblas->beta[1] = beta.imag();
    _apply<int>(_spblas->Y, _spblas->A, _spblas->X, _spblas->alpha, _spblas->beta, _spblas->common);
}

template <template <typename, typename, typename> class Matrix, typename T, typename I>
void SpBLAS<Matrix, T, I>::apply(Vector_Dense<T, I> &y, const T &alpha, const T &beta, const Vector_Dense<T, I> &x)
{
    eslog::error("SpBLAS wrapper is incompatible with T=%dB, I=%dB\n", sizeof(T), sizeof(I));
}

template <>
void SpBLAS<Matrix_CSR, float, int>::apply(Matrix_Dense<float> &y, const float &alpha, const float &beta, const Matrix_Dense<float> &x, bool trans)
{
    if (_spblas->A->nrow == 0 || _spblas->A->ncol == 0) return;
    if (trans == false) eslog::error("implement SpBLAS::apply with non-transposed matrix\n");

    update(_spblas->X, x);
    update(_spblas->Y, y);
    _spblas->alpha[0] = alpha;
    _spblas->beta[0] = beta;
    _apply<int>(_spblas->Y, _spblas->A, _spblas->X, _spblas->alpha, _spblas->beta, _spblas->common);
}

template <>
void SpBLAS<Matrix_CSR, double, int>::apply(Matrix_Dense<double> &y, const double &alpha, const double &beta, const Matrix_Dense<double> &x, bool trans)
{
    if (_spblas->A->nrow == 0 || _spblas->A->ncol == 0) return;
    if (trans == false) eslog::error("implement SpBLAS::apply with non-transposed matrix\n");

    update(_spblas->X, x);
    update(_spblas->Y, y);
    _spblas->alpha[0] = alpha;
    _spblas->beta[0] = beta;
    _apply<int>(_spblas->Y, _spblas->A, _spblas->X, _spblas->alpha, _spblas->beta, _spblas->common);
}

template <>
void SpBLAS<Matrix_CSR, std::complex<double>, int>::apply(Matrix_Dense<std::complex<double>, int> &y, const std::complex<double> &alpha, const std::complex<double> &beta, const Matrix_Dense<std::complex<double> > &x, bool trans)
{
    if (_spblas->A->nrow == 0 || _spblas->A->ncol == 0) return;
    if (trans == false) eslog::error("implement SpBLAS::apply with non-transposed matrix\n");

    update(_spblas->X, x);
    update(_spblas->Y, y);
    _spblas->alpha[0] = alpha.real();
    _spblas->alpha[1] = alpha.imag();
    _spblas->beta[0] = beta.real();
    _spblas->beta[1] = beta.imag();
    _apply<int>(_spblas->Y, _spblas->A, _spblas->X, _spblas->alpha, _spblas->beta, _spblas->common);
}

template <template <typename, typename, typename> class Matrix, typename T, typename I>
void SpBLAS<Matrix, T, I>::apply(Matrix_Dense<T, I> &y, const T &alpha, const T &beta, const Matrix_Dense<T, I> &x, bool trans)
{
    eslog::error("SpBLAS wrapper is incompatible with T=%dB, I=%dB\n", sizeof(T), sizeof(I));
}

struct _handle_trsm
{
    // cholmod_common cm_common;
    // cholmod_factor A_cm;
    // cholmod_dense X_cm;
    // int sys;
};

template<typename T, typename I>
void trsm(MatrixCsxView_new<T,I> & A, MatrixDenseView_new<T> & X, MatrixDenseView_new<T> & Y, handle_trsm & handle, char stage);
{
    eslog::error("suitesparse spblas trsm not supported\n");
    // // WARNING: does not really make sense to use
    // // better: use spsolver, factorize, and solve, there it uses supernodal which is better
    // if(A.nrows != A.ncols) eslog::error("system matrix A has to be square\n");
    // if(X.nrows != Y.nrows || X.ncols != Y.ncols) eslog::error("X and Y sizes dont match\n");
    // if(X.order != Y.order) eslog::error("X and Y orders must match\n");
    // if(A.nrows != X.nrows) eslog::error("incompatible matrix sizes\n");
    // if(A.diag == 'U') eslog::error("no support for unit diag\n");
    // if(X.order != 'C') eslog::error("only colmajor X supported");

    // if(stage == 'P') { // Preprocess
    //     _start<I>(handle->cm_common);
    //     handle->cm_common.final_ll = 1;
    //     handle->cm_common.nthreads_max = 1;
    //     handle->cm_common.nmethods = 1;
    //     handle->cm_common.method[0].ordering = CHOLMOD_METIS;
    //     handle->cm_common.itype = _getCholmodItype<I>();
    //     handle->cm_common.supernodal = CHOLMOD_SIMPLICIAL;

    //     size_t n = A.nrows;
    //     handle->A_cm.n = n;
    //     handle->A_cm.minor = n;
    //     handle->A_cm.Perm = nullptr;
    //     handle->A_cm.ColCount = nullptr;
    //     handle->A_cm.Iperm = nullptr;
    //     handle->A_cm.nzmax = A.nnz;
    //     handle->A_cm.p = A.ptrs;
    //     handle->A_cm.i = A.idxs;
    //     handle->A_cm.x = A.vals;
    //     handle->A_cm.z = nullptr;
    //     handle->A_cm.nz = malloc(n * sizeof(I));
    //     handle->A_cm.next = malloc((n + 2) * sizeof(I));
    //     handle->A_cm.prev = malloc((n + 2) * sizeof(I));
    //     handle->A_cm.ordering = -1;
    //     handle->A_cm.is_ll = 1;
    //     handle->A_cm.is_super = 0;
    //     handle->A_cm.is_monotonic = 1;
    //     handle->A_cm.itype = _getCholmodItype<I>();
    //     handle->A_cm.xtype = _getCholmodXtype<T>();
    //     handle->A_cm.dtype = _getCholmodXtype<T>();
    //     handle->A_cm.useGPU = 0;

    //     I * nz = reinterpret_cast<I*>(handle->A_cm.nz);
    //     I * next = reinterpret_cast<I*>(handle->A_cm.next);
    //     I * prev = reinterpret_cast<I*>(handle->A_cm.prev);
    //     for(size_t c = 0; c < n; c++) {
    //         nz[c] = A.ptrs[c+1] - A.ptrs[c];
    //     }
    //     for(size_t i = 0; i < n; i++) {
    //         next[i] = i+1;
    //     }
    //     next[n] = -1;
    //     next[n+1] = 0;
    //     prev[0] = n+1;
    //     for(size_t i = 1; i < n; i++) {
    //         next[i] = i-1;
    //     }
    //     prev[n] = n-1;
    //     prev[n+1] = -1;

    //     handle->X_cm.nrow = X.nrows;
    //     handle->X_cm.ncol = X.ncols;
    //     handle->X_cm.nzmax = X.get_num_blocks() * X.get_block_size();
    //     handle->X_cm.d = X.ld;
    //     handle->X_cm.x = X.vals;
    //     handle->X_cm.z = nullptr;
    //     handle->X_cm.xtype = _getCholmodXtype<T>();
    //     handle->X_cm.dtype = _getCholmodDtype<D>();

    //     handle->sys = (((A.order == 'C') == (A.uplo == 'L')) ? CHOLMOD_L : CHOLMOD_Lt);
    // }
    // if(stage == 'C') { // Compute
    //     handle->X_cm.x = X.vals;

    //     cholmod_dense * Y_cm = _solve<I>(handle->sys, handle->A_cm, handle->X_cm, handle->cm_common);
    //     MatrixDenseView_new<T> Y_tmp;
    //     Y_tmp.set(Y_cm.nrow, Y_cm.ncol, Y_cm.d, 'C', reinterpret_cast<T*>(Y_cm.x));
    //     copy_dense<T>::do_all(&Y_tmp, &Y);

    //     _free<I>(Y_cm);
    // }
    // if(stage == 'F') { // Finalize, Free
    //     free(handle->A_cm.nz);
    //     free(handle->A_cm.next);
    //     free(handle->A_cm.prev);
    //     _finish<I>(ext->cholmod.cm_common);
    // }
}

struct _handle_mm
{
    // cholmod_common cm_common;
    // cholmod_sparse A_cm;
    // cholmod_dense B_cm;
    // cholmod_dense C_cm;
};

template<typename T, typename I>
void mm(MatrixCsxView_new<T,I> & A, MatrixDenseView_new<T> & B, MatrixDenseView_new<T> & C, T alpha, T beta, handle_mm & handle, char stage)
{
    eslog::error("suitesparse spblas mm not supported\n");
    // if(A.nrows != C.nrows || B.ncols != C.ncols || A.ncols != B.nrows) eslog::error("incompatible matrices\n");
    // if(B.order != C.order) eslog::error("B and C order must match\n");

    // if(stage == 'P') {
    //     _start<I>(handle->cm_common);
    //     handle->cm_common.nthreads_max = 1;
    //     handle->cm_common.itype = _getCholmodItype<I>();



    //     handle->B_cm.nrow = B.nrows;
    //     handle->B_cm.ncol = B.ncols;
    //     handle->B_cm
        
    // }
    // if(stage == 'C') {
    //     double cm_alpha[2];
    //     double cm_beta[2];
    //     if constexpr(utils::is_real<T>()) {
    //         cm_alpha[0] = alpha;
    //         cm_alpha[1] = T{0};
    //         cm_beta[0] = beta;
    //         cm_beta[1] = T{0};
    //     }
    //     if constexpr(utils::is_complex<T>()) {
    //         cm_alpha[0] = alpha.real();
    //         cm_alpha[1] = alpha.imag();
    //         cm_beta[0] = beta.real();
    //         cm_beta[1] = beta.imag();
    //     }
    //     _apply<I>(handle->C_cm, handle->A_cm, handle->B_cm, cm_alpha, cm_beta, handle->cm_common);
    // }
    // if(stage == 'F') {
    //     _finish<I>(ext->cholmod.cm_common);
    // }
}

}

#include "math/wrappers/math.spblas.inst.hpp"

#endif
#endif
