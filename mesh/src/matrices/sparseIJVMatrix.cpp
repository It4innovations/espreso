#include "sparseIJVMatrix.h"

SparseIJVMatrix::SparseIJVMatrix(const DenseMatrix &other): Matrix(other.rows(), other.columns())
{
	MKL_INT nnz = other.nonZeroValues();
	_rowIndices.reserve(nnz);
	_columnIndices.reserve(nnz);
	_values.reserve(nnz);

	for(size_t r = 0; r < other.rows(); r++) {
		for(size_t c = 0; c < other.columns(); c++) {
			if (other(r, c) != 0) {
				_rowIndices.push_back(r);
				_columnIndices.push_back(c);
				_values.push_back(other(r, c));
			}
		}
	}
}

SparseIJVMatrix::SparseIJVMatrix(const SparseDOKMatrix &other): Matrix(other.rows(), other.columns())
{
	MKL_INT nnz = other.nonZeroValues();
	_rowIndices.reserve(nnz);
	_columnIndices.reserve(nnz);
	_values.reserve(nnz);

	const MatrixMap &dokValues = other.values();
	MatrixMap::const_iterator row;
	for(row = dokValues.begin(); row != dokValues.end(); ++row) {
		const ColumnMap &columns = row->second;

		ColumnMap::const_iterator column;
		for(column = columns.begin(); column != columns.end(); ++column) {
			if (column->second != 0) {
				_rowIndices.push_back(row->first);
				_columnIndices.push_back(column->first);
				_values.push_back(column->second);
			}
		}
	}
}

SparseIJVMatrix::SparseIJVMatrix(const SparseCSRMatrix &other): Matrix(other.rows(), other.columns())
{
	MKL_INT nnz = other.nonZeroValues();
	MKL_INT rows = _rows;
	_rowIndices.resize(nnz);
	_columnIndices.resize(nnz);
	_values.resize(nnz);

	std::cout << "NONZERO: " << nnz << "\n";

	MKL_INT job[6] = {
		0, 			// CSR to IJV
		0,			// zero based indexing
		0,			// zero based indexing
		0,			// without any meaning
		nnz,		// non-zero values
		3,			// fill all output arrays
	};

	MKL_INT info;

	mkl_dcsrcoo(
		job, &rows,
		const_cast<double*>(other.values()), const_cast<MKL_INT*>(other.columnIndices()), const_cast<MKL_INT*>(other.rowPtrs()),
		&nnz, values(), rowIndices(), columnIndices(),
		&info);
}

SparseIJVMatrix::SparseIJVMatrix(SparseVVPMatrix &other): Matrix(other.rows(), other.columns())
{
	other.shrink();
	MKL_INT nnz = other.nonZeroValues();
	_rowIndices.reserve(nnz);
	_columnIndices.reserve(nnz);
	_values.reserve(nnz);

	const VVP &values = other.values();

	for (size_t row = 0; row < _rows; row++) {
		for (size_t column = 0; column < values[row].size(); column++) {
			_rowIndices.push_back(row);
			_columnIndices.push_back(values[row][column].first);
			_values.push_back(values[row][column].second);
		}
	}
}

SparseIJVMatrix& SparseIJVMatrix::operator=(const DenseMatrix &other)
{
	SparseIJVMatrix tmp(other);
	assign(*this, tmp);
	return *this;
}

SparseIJVMatrix& SparseIJVMatrix::operator=(const SparseDOKMatrix &other)
{
	SparseIJVMatrix tmp(other);
	assign(*this, tmp);
	return *this;
}

SparseIJVMatrix& SparseIJVMatrix::operator=(const SparseCSRMatrix &other)
{
	SparseIJVMatrix tmp(other);
	assign(*this, tmp);
	return *this;
}

SparseIJVMatrix& SparseIJVMatrix::operator=(SparseVVPMatrix &other)
{
	SparseIJVMatrix tmp(other);
	assign(*this, tmp);
	return *this;
}

void SparseIJVMatrix::reserve(size_t size)
{
	_rowIndices.reserve(size);
	_columnIndices.reserve(size);
	_values.reserve(size);
}

void SparseIJVMatrix::transpose()
{
	_rowIndices.swap(_columnIndices);
}
