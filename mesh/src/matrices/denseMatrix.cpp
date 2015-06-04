#include "denseMatrix.h"

DenseMatrix::DenseMatrix(const SparseDOKMatrix &other): Matrix(other.rows(), other.columns(), DenseMatrixIndexing)
{
	std::vector<double>(rows() * columns(), 0).swap(_values);

	MatrixMap::const_iterator row;
	for(row = other.values().begin(); row != other.values().end(); ++row) {
		ColumnMap::const_iterator column;
		for(column = row->second.begin(); column != row->second.end(); ++column) {
			set(row->first - other.indexing(), column->first - other.indexing(), column->second);
		}
	}
}

DenseMatrix::DenseMatrix(const SparseCSRMatrix &other): Matrix(other.rows(), other.columns(), DenseMatrixIndexing)
{
	std::vector<double>(rows() * columns(), 0).swap(_values);

	const eslocal *rowPtrs = other.rowPtrs();
	const eslocal *columnIndices = other.columnIndices();
	const double *values = other.values();

	for (eslocal r = 0; r < _rows; r++) {
		for (eslocal i = rowPtrs[r]; i < rowPtrs[r + 1]; i++) {
			set(r, columnIndices[i - other.indexing()] - other.indexing(), values[i - other.indexing()]);
		}
	}
}

DenseMatrix::DenseMatrix(const SparseIJVMatrix &other): Matrix(other.rows(), other.columns(), DenseMatrixIndexing)
{
	std::vector<double>(rows() * columns(), 0).swap(_values);

	const eslocal *rowIndices = other.rowIndices();
	const eslocal *columnIndices = other.columnIndices();
	const double *values = other.values();

	for (eslocal i = 0; i < other.nonZeroValues(); i++) {
		set(rowIndices[i] - other.indexing(), columnIndices[i] - other.indexing(), values[i]);
	}
}

DenseMatrix& DenseMatrix::operator=(const SparseDOKMatrix &other)
{
	DenseMatrix tmp(other);
	assign(*this, tmp);
	return *this;
}

DenseMatrix& DenseMatrix::operator=(const SparseCSRMatrix &other)
{
	DenseMatrix tmp(other);
	assign(*this, tmp);
	return *this;
}

DenseMatrix& DenseMatrix::operator=(const SparseIJVMatrix &other)
{
	DenseMatrix tmp(other);
	assign(*this, tmp);
	return *this;
}

DenseMatrix& DenseMatrix::operator=(double value)
{
	std::fill(_values.begin(), _values.end(), value);
	return *this;
}

DenseMatrix DenseMatrix::operator*(DenseMatrix &M)
{
	DenseMatrix result;
	result.multiply(*this, M);
	return result;
}

void DenseMatrix::multiply(
		const DenseMatrix &A, const DenseMatrix &B,
		double alfa, double beta,
		bool transposeA, bool transposeB)
{
	if ( (transposeA ? A.rows() : A.columns()) != (transposeB ? B.columns() : B.rows()) ) {
		std::cerr << "Matrix multiplication: matrices have incorrect dimensions.\n";
		exit(EXIT_FAILURE);
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


