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

	const MKL_INT *rowPtrs = other.rowPtrs();
	const MKL_INT *columnIndices = other.columnIndices();
	const double *values = other.values();

	for (MKL_INT r = 0; r < _rows; r++) {
		for (MKL_INT i = rowPtrs[r]; i < rowPtrs[r + 1]; i++) {
			set(r, columnIndices[i - other.indexing()] - other.indexing(), values[i - other.indexing()]);
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
