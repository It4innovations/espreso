
#include "math/wrappers/math.spblas.h"
#include "esinfo/eslog.h"

#include <complex>

#ifndef HAVE_MKL
#ifdef HAVE_SUITESPARSE

#include "wrappers/suitesparse/w.suitesparse.cholmod.h"

namespace espreso {

struct Matrix_SpBLAS_External_Representation {
	cholmod_sparse *A;
	cholmod_dense *X, *Y;
	double alpha[2] = { 1., 1. }, beta[2] = { 0., 0. };
	cholmod_common common;

	Matrix_SpBLAS_External_Representation(): A(new cholmod_sparse()), X(new cholmod_dense()), Y(new cholmod_dense()) { }
	~Matrix_SpBLAS_External_Representation() { delete A; delete X; delete Y;}
};

template <typename T, template <typename> class Matrix>
SpBLAS<T, Matrix>::SpBLAS()
: matrix{}, _spblas{}
{

}

template <typename T, template <typename> class Matrix>
SpBLAS<T, Matrix>::~SpBLAS()
{
	if (_spblas) {
		_finish<esint>(_spblas->common);
		delete _spblas;
	}
}

template <typename T, template <typename> class Matrix>
SpBLAS<T, Matrix>::SpBLAS(Matrix<T> &a)
: matrix(&a)
{
	_spblas = new Matrix_SpBLAS_External_Representation();
	_start<esint>(_spblas->common);
	setSymmetric(_spblas->A, *matrix);
	updateSymmetric(_spblas->A, *matrix);
}

template <typename T, template <typename> class Matrix>
void SpBLAS<T, Matrix>::insert(Matrix<T> &a)
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

template <typename T, template <typename> class Matrix>
void SpBLAS<T, Matrix>::insertTransposed(Matrix<T> &a)
{
	if (_spblas) {
		_finish<esint>(_spblas->common);
		delete _spblas;
	}
	_spblas = new Matrix_SpBLAS_External_Representation();
	_start<esint>(_spblas->common);
	setSymmetric(_spblas->A, a);
	updateSymmetric(_spblas->A, a);
}

template <typename T, template <typename> class Matrix>
void SpBLAS<T, Matrix>::extractUpper(Matrix<T> &a)
{
	_extractUpper(_spblas->A, _spblas->common, a);
}

template <>
void SpBLAS<float, Matrix_CSR>::apply(Vector_Dense<float> &y, const float &alpha, const float &beta, const Vector_Dense<float> &x)
{
	update(_spblas->X, x);
	update(_spblas->Y, y);
	_spblas->alpha[0] = alpha;
	_spblas->beta[0] = beta;
	_apply<esint>(_spblas->Y, _spblas->A, _spblas->X, _spblas->alpha, _spblas->beta, _spblas->common);
}

template <>
void SpBLAS<double, Matrix_CSR>::apply(Vector_Dense<double> &y, const double &alpha, const double &beta, const Vector_Dense<double> &x)
{
	update(_spblas->X, x);
	update(_spblas->Y, y);
	_spblas->alpha[0] = alpha;
	_spblas->beta[0] = beta;
	_apply<esint>(_spblas->Y, _spblas->A, _spblas->X, _spblas->alpha, _spblas->beta, _spblas->common);
}

template <>
void SpBLAS<std::complex<double>, Matrix_CSR>::apply(Vector_Dense<std::complex<double>> &y, const std::complex<double> &alpha, const std::complex<double> &beta, const Vector_Dense<std::complex<double>> &x)
{
	update(_spblas->X, x);
	update(_spblas->Y, y);
	_spblas->alpha[0] = alpha.real();
	_spblas->alpha[1] = alpha.imag();
	_spblas->beta[0] = beta.real();
	_spblas->beta[1] = beta.imag();
	_apply<esint>(_spblas->Y, _spblas->A, _spblas->X, _spblas->alpha, _spblas->beta, _spblas->common);
}

template <typename T, template <typename> class Matrix>
void SpBLAS<T, Matrix>::transposeTo(SpBLAS<T, Matrix> &A)
{
	if (A._spblas) {
		_finish<esint>(A._spblas->common);
		delete _spblas;
	}
	A._spblas = new Matrix_SpBLAS_External_Representation();
	_start<esint>(A._spblas->common);
	_transpose<esint>(_spblas->A, A._spblas->A, _spblas->common);
	A.matrix = nullptr;
}

template <>
void SpBLAS<float, Matrix_CSR>::multiply(SpBLAS<float, Matrix_CSR> &A, SpBLAS<float, Matrix_CSR> &B)
{
	if (_spblas) {
		_finish<esint>(_spblas->common);
		delete _spblas;
	}
	_spblas = new Matrix_SpBLAS_External_Representation();
	_start<esint>(_spblas->common);
	_multiply<esint>(A._spblas->A, B._spblas->A, _spblas->A, _spblas->common);
}

template <>
void SpBLAS<double, Matrix_CSR>::multiply(SpBLAS<double, Matrix_CSR> &A, SpBLAS<double, Matrix_CSR> &B)
{
	if (_spblas) {
		_finish<esint>(_spblas->common);
		delete _spblas;
	}
	_spblas = new Matrix_SpBLAS_External_Representation();
	_start<esint>(_spblas->common);
	_multiply<esint>(A._spblas->A, B._spblas->A, _spblas->A, _spblas->common);
}

template struct SpBLAS<float, Matrix_CSR>;
template struct SpBLAS<double, Matrix_CSR>;
template struct SpBLAS<std::complex<double>, Matrix_CSR>;

}

#endif
#endif
