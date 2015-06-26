#include "sparseIJVMatrix.h"

std::ostream& operator<<(std::ostream& os, const SparseIJVMatrix &m)
{
  os << m._rows << " " << m._columns << " " << m.nonZeroValues() << "\n";

	for (size_t i = 0; i < m.nonZeroValues(); i++) {
		os << m._rowIndices[i] << " " << m._columnIndices[i] << " " << m._values[i] << "\n";
	}
	return os;
}

SparseIJVMatrix::SparseIJVMatrix(const DenseMatrix &other): Matrix(other.rows(), other.columns(), IJVMatrixIndexing)
{
	eslocal nnz = other.nonZeroValues();
	_rowIndices.reserve(nnz);
	_columnIndices.reserve(nnz);
	_values.reserve(nnz);

	for(size_t r = 0; r < other.rows(); r++) {
		for(size_t c = 0; c < other.columns(); c++) {
			if (Matrix::nonZero(other(r, c))) {
				_rowIndices.push_back(r + _indexing);
				_columnIndices.push_back(c + _indexing);
				_values.push_back(other(r, c));
			}
		}
	}
}

SparseIJVMatrix::SparseIJVMatrix(const SparseDOKMatrix &other): Matrix(other.rows(), other.columns(), IJVMatrixIndexing)
{
	eslocal nnz = other.nonZeroValues();
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
				_rowIndices.push_back(row->first - other.indexing() + _indexing);
				_columnIndices.push_back(column->first - other.indexing() + _indexing);
				_values.push_back(column->second);
			}
		}
	}
}

SparseIJVMatrix::SparseIJVMatrix(const SparseCSRMatrix &other): Matrix(other.rows(), other.columns(), IJVMatrixIndexing)
{
	eslocal nnz = other.nonZeroValues();
	eslocal rows = _rows;
	_rowIndices.resize(nnz);
	_columnIndices.resize(nnz);
	_values.resize(nnz);

	eslocal job[6] = {
		0, 					// CSR to IJV
		other.indexing(),	// indexing of CSR matrix
		_indexing,			// indexing of IJV matrix
		0,					// without any meaning
		nnz,				// non-zero values
		3,					// fill all output arrays
	};

	eslocal info;

	mkl_dcsrcoo(
		job, &rows,
		const_cast<double*>(other.values()), const_cast<eslocal*>(other.columnIndices()), const_cast<eslocal*>(other.rowPtrs()),
		&nnz, values(), rowIndices(), columnIndices(),
		&info);
}

SparseIJVMatrix::SparseIJVMatrix(SparseVVPMatrix &other): Matrix(other.rows(), other.columns(), IJVMatrixIndexing)
{
	other.shrink();
	eslocal nnz = other.nonZeroValues();
	_rowIndices.reserve(nnz);
	_columnIndices.reserve(nnz);
	_values.reserve(nnz);

	const VVP &values = other.values();

	for (size_t row = 0; row < other.rows(); row++) {
		for (size_t column = 0; column < values[row].size(); column++) {
			_rowIndices.push_back(row + _indexing);
			_columnIndices.push_back(values[row][column].first - other.indexing() + _indexing);
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
