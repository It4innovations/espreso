
#include "math/wrappers/math.spblas.h"
#include "esinfo/eslog.h"

#include <complex>

#ifndef HAVE_MKL
#ifdef HAVE_SUITESPARSE

#include "wrappers/suitesparse/w.suitesparse.cholmod.h"

namespace espreso {

template struct SpBLAS<Matrix_CSR, float, int>;
template struct SpBLAS<Matrix_CSR, double, int>;
template struct SpBLAS<Matrix_CSR, std::complex<double>, int>;

struct Matrix_SpBLAS_External_Representation {
	cholmod_sparse *A;
	cholmod_dense *X, *Y;
	double alpha[2] = { 1., 1. }, beta[2] = { 0., 0. };
	cholmod_common common;

	Matrix_SpBLAS_External_Representation(): A(new cholmod_sparse()), X(new cholmod_dense()), Y(new cholmod_dense()) { }
	~Matrix_SpBLAS_External_Representation() { delete A; delete X; delete Y;}
};

template <template <typename, typename> class Matrix, typename T, typename I>
SpBLAS<Matrix, T, I>::SpBLAS()
: matrix{}, _spblas{}
{

}

template <template <typename, typename> class Matrix, typename T, typename I>
SpBLAS<Matrix, T, I>::~SpBLAS()
{
	if (_spblas) {
		_finish<esint>(_spblas->common);
		delete _spblas;
	}
}

template <template <typename, typename> class Matrix, typename T, typename I>
SpBLAS<Matrix, T, I>::SpBLAS(Matrix<T, I> &a)
: matrix(&a)
{
	_spblas = new Matrix_SpBLAS_External_Representation();
	_start<esint>(_spblas->common);
	setSymmetric(_spblas->A, *matrix);
	updateSymmetric(_spblas->A, *matrix);
}

template <template <typename, typename> class Matrix, typename T, typename I>
void SpBLAS<Matrix, T, I>::insert(Matrix<T, I> &a)
{
	matrix = &a;
	if (_spblas) {
		_finish<esint>(_spblas->common);
		delete _spblas;
	}
	_spblas = new Matrix_SpBLAS_External_Representation();
	_start<esint>(_spblas->common);
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
	_apply<esint>(_spblas->Y, _spblas->A, _spblas->X, _spblas->alpha, _spblas->beta, _spblas->common);
}

template <>
void SpBLAS<Matrix_CSR, double, int>::apply(Vector_Dense<double> &y, const double &alpha, const double &beta, const Vector_Dense<double> &x)
{
	update(_spblas->X, x);
	update(_spblas->Y, y);
	_spblas->alpha[0] = alpha;
	_spblas->beta[0] = beta;
	_apply<esint>(_spblas->Y, _spblas->A, _spblas->X, _spblas->alpha, _spblas->beta, _spblas->common);
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
	_apply<esint>(_spblas->Y, _spblas->A, _spblas->X, _spblas->alpha, _spblas->beta, _spblas->common);
}

}

#endif
#endif
