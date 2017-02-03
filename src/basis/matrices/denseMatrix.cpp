#include "denseMatrix.h"

using namespace espreso;

DenseMatrix& DenseMatrix::operator=(double value)
{
	std::fill(_values.begin(), _values.end(), value);
	return *this;
}

DenseMatrix& DenseMatrix::operator+=(const DenseMatrix &other)
{
	for (size_t r = 0; r < other.rows(); r++) {
		for (size_t c = 0; c < other.columns(); c++) {
			(*this)(r, c) += other(r, c);
		}
	}
	return *this;
}

DenseMatrix DenseMatrix::operator*(DenseMatrix &other)
{
	DenseMatrix result;
	result.multiply(*this, other);
	return result;
}

void DenseMatrix::multiply(
		const DenseMatrix &A, const DenseMatrix &B,
		double alfa, double beta,
		bool transposeA, bool transposeB)
{
	if ( (transposeA ? A.rows() : A.columns()) != (transposeB ? B.columns() : B.rows()) ) {
		ESINFO(ERROR) << "Matrix multiplication: matrices have incorrect dimensions.";
	}
	resize(transposeA ? A.columns() : A.rows(), transposeB ? B.rows() : B.columns());

	cblas_dgemm(
		CblasRowMajor,
		transposeA ? CblasTrans : CblasNoTrans,
		transposeB ? CblasTrans : CblasNoTrans,
		transposeA ? A.columns() : A.rows(),
		transposeB ? B.rows() : B.columns(),
		transposeA ? A.rows() : A.columns(),
		alfa,
		A.values(), A.columns(),
		B.values(), B.columns(),
		beta,
		values(), transposeB ? B.rows() : B.columns());
}

void DenseMatrix::resize(size_t rows, size_t columns)
{
	_rows = rows;
	_columns = columns;
	_values.resize(rows * columns);
}

void DenseMatrix::transpose()
{
	std::vector<double> copy(_values.size());

	MKL_Domatcopy(
			'r', 't',
			_rows, _columns,
			1, values(), _columns,
			&copy[0], _rows);

	_values.swap(copy);
	size_t tmp = _rows;
	_rows = _columns;
	_columns = tmp;
}


