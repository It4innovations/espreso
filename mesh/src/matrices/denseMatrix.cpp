#include "denseMatrix.h"

DenseMatrix::DenseMatrix(const SparseDOKMatrix &other): Matrix(other.rows(), other.columns(), DenseMatrixIndexing)
{
	std::vector<double>(rows() * columns(), 0).swap(_values);

	MatrixMap::const_iterator row;
	for(row = other.values().begin(); row != other.values().end(); ++row) {
		ColumnMap::const_iterator column;
		for(column = row->second.begin(); column != row->second.end(); ++column) {
			this->operator ()(row->first, column->first) = column->second;
		}
	}
}

DenseMatrix::DenseMatrix(const SparseCSRMatrix &other): Matrix(other.rows(), other.columns(), DenseMatrixIndexing)
{
	std::vector<double>(rows() * columns(), 0).swap(_values);

	const MKL_INT *rowPtrs = other.rowPtrs();
	const MKL_INT *columnIndices = other.columnIndices();
	const double *values = other.values();

	for (MKL_INT r = 0; r < _rows; r++) {
		for (MKL_INT i = rowPtrs[r]; i < rowPtrs[r + 1]; i++) {
			this->operator ()(r, columnIndices[i]) = values[i];
		}
	}
}

DenseMatrix::DenseMatrix(const SparseIJVMatrix &other): Matrix(other.rows(), other.columns(), DenseMatrixIndexing)
{
	std::vector<double>(rows() * columns(), 0).swap(_values);

	const MKL_INT *rowIndices = other.rowIndices();
	const MKL_INT *columnIndices = other.columnIndices();
	const double *values = other.values();

	for (MKL_INT i = 0; i < other.nonZeroValues(); i++) {
		this->operator ()(rowIndices[i], columnIndices[i]) = values[i];
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

void DenseMatrix::resize(size_t rows, size_t columns)
{
	_rows = rows;
	_columns = columns;
	_values.resize(rows * columns);
}

void DenseMatrix::transpose()
{
	std::vector<double> values(_rows * _columns);
	for (MKL_INT r = 0; r < _rows; r++) {
		for (MKL_INT c = 0; c < _columns; c++) {
			values[c * _rows + r] = _values[r * _columns + c];
		}
	}

	_values.swap(values);
}
