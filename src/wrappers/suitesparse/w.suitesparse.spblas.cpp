
#include "math/wrappers/math.spblas.h"
#include "esinfo/eslog.h"

#include <complex>

#ifndef HAVE_MKL
#ifdef HAVE_SUITESPARSE

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
SpBLAS<Matrix, T, I>::SpBLAS(MatrixType &a)
: matrix(&a)
{
    _spblas = new Matrix_SpBLAS_External_Representation();
    _start<I>(_spblas->common);
    setSymmetric(_spblas->A, *matrix);
    updateSymmetric(_spblas->A, *matrix);
}

template <template <typename, typename, typename> class Matrix, typename T, typename I>
void SpBLAS<Matrix, T, I>::insert(MatrixType &a)
{
    matrix = &a;
    if (_spblas) {
        _finish<I>(_spblas->common);
        delete _spblas;
    }
    _spblas = new Matrix_SpBLAS_External_Representation();
    _start<I>(_spblas->common);
    if (a.shape == Matrix_Shape::FULL) {
        setAsymmetric(_spblas->A, _spblas->common, *matrix);
        matrix = nullptr;
    } else {
        setSymmetric(_spblas->A, *matrix);
        updateSymmetric(_spblas->A, *matrix);
    }
}

template <>
void SpBLAS<Matrix_CSR, float, int>::apply(Vector_Dense<float> &y, const float &alpha, const float &beta, const Vector_Dense<float> &x)
{
    update(_spblas->X, x);
    update(_spblas->Y, y);
    _spblas->alpha[0] = alpha;
    _spblas->beta[0] = beta;
    _apply<int>(_spblas->Y, _spblas->A, _spblas->X, _spblas->alpha, _spblas->beta, _spblas->common);
}

template <>
void SpBLAS<Matrix_CSR, double, int>::apply(Vector_Dense<double> &y, const double &alpha, const double &beta, const Vector_Dense<double> &x)
{
    update(_spblas->X, x);
    update(_spblas->Y, y);
    _spblas->alpha[0] = alpha;
    _spblas->beta[0] = beta;
    _apply<int>(_spblas->Y, _spblas->A, _spblas->X, _spblas->alpha, _spblas->beta, _spblas->common);
}

template <>
void SpBLAS<Matrix_CSR, std::complex<double>, int>::apply(Vector_Dense<std::complex<double>, int> &y, const std::complex<double> &alpha, const std::complex<double> &beta, const Vector_Dense<std::complex<double> > &x)
{
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

}

#include "math/wrappers/math.spblas_inst.hpp"

#endif
#endif
