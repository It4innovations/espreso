
#include "math/wrappers/math.spblas.h"
#include "esinfo/eslog.h"

#include <complex>

#ifdef HAVE_SUITESPARSE
#ifdef ESPRESO_USE_WRAPPER_SPBLAS_SUITESPARSE

#include "w.suitesparse.cholmod.h"

#include "math/operations/convert_dnx_dny.h"

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
    _apply<int>(_spblas->Y, _spblas->A, _spblas->X, false, _spblas->alpha, _spblas->beta, _spblas->common);
}

template <>
void SpBLAS<Matrix_CSR, double, int>::apply(Vector_Dense<double> &y, const double &alpha, const double &beta, const Vector_Dense<double> &x)
{
    if (_spblas->A->nrow == 0 || _spblas->A->ncol == 0) return;

    update(_spblas->X, x);
    update(_spblas->Y, y);
    _spblas->alpha[0] = alpha;
    _spblas->beta[0] = beta;
    _apply<int>(_spblas->Y, _spblas->A, _spblas->X, false, _spblas->alpha, _spblas->beta, _spblas->common);
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
    _apply<int>(_spblas->Y, _spblas->A, _spblas->X, false, _spblas->alpha, _spblas->beta, _spblas->common);
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
    _apply<int>(_spblas->Y, _spblas->A, _spblas->X, false, _spblas->alpha, _spblas->beta, _spblas->common);
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
    _apply<int>(_spblas->Y, _spblas->A, _spblas->X, false, _spblas->alpha, _spblas->beta, _spblas->common);
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
    _apply<int>(_spblas->Y, _spblas->A, _spblas->X, false, _spblas->alpha, _spblas->beta, _spblas->common);
}

template <template <typename, typename, typename> class Matrix, typename T, typename I>
void SpBLAS<Matrix, T, I>::apply(Matrix_Dense<T, I> &y, const T &alpha, const T &beta, const Matrix_Dense<T, I> &x, bool trans)
{
    eslog::error("SpBLAS wrapper is incompatible with T=%dB, I=%dB\n", sizeof(T), sizeof(I));
}

namespace math {
namespace spblas {

struct _handle_trsm
{
};

template<typename T, typename I>
void trsm(MatrixCsxView_new<T,I> & A, MatrixDenseView_new<T> & X, MatrixDenseView_new<T> & Y, handle_trsm & handle, char stage)
{
    eslog::error("suitesparse spblas trsm not supported yet\n");
}

struct _handle_mm
{
};

template<typename T, typename I>
void mm(MatrixCsxView_new<T,I> & A, MatrixDenseView_new<T> & B, MatrixDenseView_new<T> & C, T alpha, T beta, handle_mm & handle, char stage)
{
    if(A.nrows != C.nrows || B.ncols != C.ncols || A.ncols != B.nrows) eslog::error("incompatible matrices\n");

    if(A.nrows == 0) return;
    if(B.ncols == 0) return;

    if(stage == 'A') { // All at once
        // mm<T,I>(A, B, C, alpha, beta, handle, 'P');
        mm<T,I>(A, B, C, alpha, beta, handle, 'C');
    }
    if(stage == 'P') { // Preprocess
    }
    if(stage == 'C') { // Compute
        if(B.order != 'C') {
            MatrixDenseData_new<T> B2;
            B2.set(B.nrows, B.ncols, 'C', AllocatorCPU_new::get_singleton());
            B2.alloc();
            math::operations::convert_dnx_dny<T>::do_all(&B, &B2, false);
            mm<T,I>(A, B2, C, alpha, beta, handle, stage);
            return;
        }
        if(C.order != 'C') {
            MatrixDenseData_new<T> C2;
            C2.set(C.nrows, C.ncols, 'C', AllocatorCPU_new::get_singleton());
            C2.alloc();
            if(beta != T{0}) math::operations::convert_dnx_dny<T>::do_all(&C, &C2, false);
            mm<T,I>(A, B, C2, alpha, beta, handle, stage);
            math::operations::convert_dnx_dny<T>::do_all(&C2, &C, false);
            return;
        }
        bool trans_A = (A.order == 'R');
        MatrixCsxView_new<T,I> A_rt = A.get_transposed_reordered_view();

        cholmod_sparse cm_A;
        cholmod_dense cm_B;
        cholmod_dense cm_C;
        if(A.order == 'R') popullateCholmodSparse<T,I>(cm_A, A_rt);
        if(A.order == 'C') popullateCholmodSparse<T,I>(cm_A, A);
        popullateCholmodDense<T,I>(cm_B, B);
        popullateCholmodDense<T,I>(cm_C, C);

        double alpha_[2];
        double beta_[2];
        _getCholmodScalar(alpha, alpha_);
        _getCholmodScalar(beta, beta_);

        cholmod_common cm_common;
        _start<I>(cm_common);
        cm_common.nthreads_max = 1;
        cm_common.itype = _getCholmodItype<I>();

        _apply<I>(&cm_C, &cm_A, &cm_B, trans_A, alpha_, beta_, cm_common);

        _finish<I>(cm_common);
    }
}

}
}

}

#include "math/wrappers/math.spblas.inst.hpp"

#endif
#endif
