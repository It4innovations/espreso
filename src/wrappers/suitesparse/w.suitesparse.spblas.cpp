
#include "math/wrappers/math.spblas.h"
#include "esinfo/eslog.h"

#include <complex>

#ifdef HAVE_SUITESPARSE
#ifdef USE_SPBLAS_SUITESPARSE

#include "wrappers/suitesparse/w.suitesparse.cholmod.h"

namespace espreso {

struct Matrix_CSR_External_Representation {
	cholmod_sparse *A;
	cholmod_dense *X, *Y;
	double alpha[2] = { 1., 1. }, beta[2] = { 0., 0. };
	cholmod_common common;

	Matrix_CSR_External_Representation(): A(new cholmod_sparse()), X(new cholmod_dense()), Y(new cholmod_dense()) { }
	~Matrix_CSR_External_Representation() { delete A; delete X; delete Y;}
};

struct Matrix_IJV_External_Representation {

	Matrix_IJV_External_Representation()
	{
		eslog::error("IJV matrix operations are not supported.\n");
	}
};

namespace math {



template <>
void commit(Matrix_Dense<double> &A)
{

}

template <>
void commit(Matrix_Dense<std::complex<double> > &A)
{

}

template <>
void commit(Matrix_CSR<double> &A)
{
	A._spblas = new Matrix_CSR_External_Representation();
	_start<esint>(A._spblas->common);
	set(A._spblas->A, A);
	update(A._spblas->A, A);
}

template <>
void commit(Matrix_CSR<std::complex<double> > &A)
{
	A._spblas = new Matrix_CSR_External_Representation();
	_start<esint>(A._spblas->common);
	set(A._spblas->A, A);
	update(A._spblas->A, A);
}

template <>
void commit(Matrix_IJV<double> &A)
{
	A._spblas = new Matrix_IJV_External_Representation();
}

template <>
void commit(Matrix_IJV<std::complex<double> > &A)
{
	A._spblas = new Matrix_IJV_External_Representation();
}

template <>
void free(Matrix_Dense<double> &A)
{

}

template <>
void free(Matrix_Dense<std::complex<double> > &A)
{

}

template <>
void free(Matrix_CSR<double> &A)
{
	if (A._spblas) {
		_finish<esint>(A._spblas->common);
		delete A._spblas;
		A._spblas = nullptr;
	}
}

template <>
void free(Matrix_CSR<std::complex<double> > &A)
{
	if (A._spblas) {
		_finish<esint>(A._spblas->common);
		delete A._spblas;
		A._spblas = nullptr;
	}
}

template <>
void free(Matrix_IJV<double> &A)
{
	if (A._spblas) {
		delete A._spblas;
		A._spblas = nullptr;
	}
}

template <>
void free(Matrix_IJV<std::complex<double> > &A)
{
	if (A._spblas) {
		delete A._spblas;
		A._spblas = nullptr;
	}
}


template <>
void apply(Vector_Dense<double> &y, const double &alpha, const Matrix_CSR<double> &A, const double &beta, const Vector_Dense<double> &x)
{
	update(A._spblas->X, x);
	update(A._spblas->Y, y);
	A._spblas->alpha[0] = alpha;
	A._spblas->beta[0] = beta;
	_apply<esint>(A._spblas->Y, A._spblas->A, A._spblas->X, A._spblas->alpha, A._spblas->beta, A._spblas->common);
}

template <>
void apply(Vector_Dense<std::complex<double> > &y, const std::complex<double> &alpha, const Matrix_CSR<std::complex<double> > &A, const std::complex<double> &beta, const Vector_Dense<std::complex<double> > &x)
{
	update(A._spblas->X, x);
	update(A._spblas->Y, y);
	A._spblas->alpha[0] = alpha.real();
	A._spblas->alpha[1] = alpha.imag();
	A._spblas->beta[0] = beta.real();
	A._spblas->beta[1] = beta.imag();
	_apply<esint>(A._spblas->Y, A._spblas->A, A._spblas->X, A._spblas->alpha, A._spblas->beta, A._spblas->common);
}

template <>
void apply(Vector_Dense<double> &y, const double &alpha, const Matrix_IJV<double> &A, const double &beta, const Vector_Dense<double> &x)
{
	eslog::error("IJV matrix operations are not supported.\n");
}

template <>
void apply(Vector_Dense<std::complex<double> > &y, const std::complex<double> &alpha, const Matrix_IJV<std::complex<double> > &A, const std::complex<double> &beta, const Vector_Dense<std::complex<double> > &x)
{
	eslog::error("IJV matrix operations are not supported.\n");
}


template <typename T> void submatrix(const Matrix_CSR<T> &input, Matrix_Dense<T> &output, esint start_row, esint end_row, esint start_col, esint end_col, bool trans, bool conj, bool output_force_full)
{
	eslog::error("Implement extract block CSR->Dense in suitesparse spblas wrapper.\n");
}

template <typename T> void submatrix(const Matrix_CSR<T> &input, Matrix_CSR<T>   &output, esint start_row, esint end_row, esint start_col, esint end_col, bool trans, bool conj, bool output_force_full)
{
	eslog::error("Implement extract block CSR->CSR in suitesparse spblas wrapper.\n");
}

template void submatrix<double>(const Matrix_CSR<double> &, Matrix_Dense<double> &, esint, esint, esint, esint, bool, bool, bool);
template void submatrix<double>(const Matrix_CSR<double> &, Matrix_CSR<double> &,   esint, esint, esint, esint, bool, bool, bool);
template void submatrix<std::complex<double>>(const Matrix_CSR<std::complex<double>> &, Matrix_Dense<std::complex<double>> &, esint, esint, esint, esint, bool, bool, bool);
template void submatrix<std::complex<double>>(const Matrix_CSR<std::complex<double>> &, Matrix_CSR<std::complex<double>> &,   esint, esint, esint, esint, bool, bool, bool);

}
}

#endif
#endif
