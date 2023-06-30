
#include "math/wrappers/math.spblas.h"
#include "esinfo/eslog.h"

#include <complex>

#ifdef HAVE_SUITESPARSE
#ifdef USE_SPBLAS_SUITESPARSE

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
SpBLAS<T, Matrix>::SpBLAS(const Matrix<T> &a)
: matrix(&a)
{
	_spblas = new Matrix_SpBLAS_External_Representation();
	_start<esint>(_spblas->common);
	set(_spblas->A, *matrix);
	update(_spblas->A, *matrix);
}

template <typename T, template <typename> class Matrix>
void SpBLAS<T, Matrix>::commit(const Matrix<T> &a)
{
	matrix = &a;
	if (_spblas) {
		_finish<esint>(_spblas->common);
		delete _spblas;
	}
	_spblas = new Matrix_SpBLAS_External_Representation();
	_start<esint>(_spblas->common);
	set(_spblas->A, *matrix);
	update(_spblas->A, *matrix);
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

template class SpBLAS<double, Matrix_CSR>;
template class SpBLAS<std::complex<double>, Matrix_CSR>;

}

#endif
#endif
