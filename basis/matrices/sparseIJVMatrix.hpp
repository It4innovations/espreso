
#include "sparseIJVMatrix.h"

namespace espreso {

template<typename Tindices>
std::ostream& operator<<(std::ostream& os, const SparseIJVMatrix<Tindices> &m)
{
	os << m.rows() << " " << m.columns() << " " << m.nonZeroValues() << "\n";

	for (size_t i = 0; i < m.nonZeroValues(); i++) {
		os << m._rowIndices[i] << " " << m._columnIndices[i] << " " << m._values[i] << "\n";
	}
	return os;
}

template<typename Tindices>
SparseIJVMatrix<Tindices>::SparseIJVMatrix(const DenseMatrix &other): Matrix(other.rows(), other.columns(), IJVMatrixIndexing)
{
	size_t nnz = other.nonZeroValues();
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

template<typename Tindices>
SparseIJVMatrix<Tindices>::SparseIJVMatrix(const SparseDOKMatrix<Tindices> &other): Matrix(other.rows(), other.columns(), IJVMatrixIndexing)
{
	size_t nnz = other.nonZeroValues();
	_rowIndices.reserve(nnz);
	_columnIndices.reserve(nnz);
	_values.reserve(nnz);

	const std::map<Tindices, std::map<Tindices, double> > &dokValues = other.values();
	typename std::map<Tindices, std::map<Tindices, double> >::const_iterator row;
	for(row = dokValues.begin(); row != dokValues.end(); ++row) {
		const std::map<Tindices, double> &columns = row->second;

		typename std::map<Tindices, double>::const_iterator column;
		for(column = columns.begin(); column != columns.end(); ++column) {
			if (column->second != 0) {
				_rowIndices.push_back(row->first - other.indexing() + _indexing);
				_columnIndices.push_back(column->first - other.indexing() + _indexing);
				_values.push_back(column->second);
			}
		}
	}
}

template<typename Tindices>
SparseIJVMatrix<Tindices>::SparseIJVMatrix(const SparseCSRMatrix<Tindices> &other): Matrix(other.rows(), other.columns(), IJVMatrixIndexing)
{
	Tindices nnz = other.nonZeroValues();
	Tindices rows = _rows;
	_rowIndices.resize(nnz);
	_columnIndices.resize(nnz);
	_values.resize(nnz);

	Tindices job[6] = {
		0, 					// CSR to IJV
		other.indexing(),	// indexing of CSR matrix
		_indexing,			// indexing of IJV matrix
		0,					// without any meaning
		nnz,				// non-zero values
		3,					// fill all output arrays
	};

	Tindices info;

	mkl_dcsrcoo(
		job, &rows,
		const_cast<double*>(other.values()), const_cast<Tindices*>(other.columnIndices()), const_cast<Tindices*>(other.rowPtrs()),
		&nnz, values(), rowIndices(), columnIndices(),
		&info);
}

template<typename Tindices>
SparseIJVMatrix<Tindices>::SparseIJVMatrix(SparseVVPMatrix<Tindices> &other): Matrix(other.rows(), other.columns(), IJVMatrixIndexing)
{
	other.shrink();
	size_t nnz = other.nonZeroValues();
	_rowIndices.reserve(nnz);
	_columnIndices.reserve(nnz);
	_values.reserve(nnz);

	const VVP<Tindices> &values = other.values();

	for (size_t row = 0; row < other.rows(); row++) {
		for (size_t column = 0; column < values[row].size(); column++) {
			_rowIndices.push_back(row + _indexing);
			_columnIndices.push_back(values[row][column].first - other.indexing() + _indexing);
			_values.push_back(values[row][column].second);
		}
	}
}

template<typename Tindices>
SparseIJVMatrix<Tindices>& SparseIJVMatrix<Tindices>::operator=(const DenseMatrix &other)
{
	SparseIJVMatrix<Tindices> tmp(other);
	assign(*this, tmp);
	return *this;
}

template<typename Tindices>
SparseIJVMatrix<Tindices>& SparseIJVMatrix<Tindices>::operator=(const SparseDOKMatrix<Tindices> &other)
{
	SparseIJVMatrix<Tindices> tmp(other);
	assign(*this, tmp);
	return *this;
}

template<typename Tindices>
SparseIJVMatrix<Tindices>& SparseIJVMatrix<Tindices>::operator=(const SparseCSRMatrix<Tindices> &other)
{
	SparseIJVMatrix<Tindices> tmp(other);
	assign(*this, tmp);
	return *this;
}

template<typename Tindices>
SparseIJVMatrix<Tindices>& SparseIJVMatrix<Tindices>::operator=(SparseVVPMatrix<Tindices> &other)
{
	SparseIJVMatrix<Tindices> tmp(other);
	assign(*this, tmp);
	return *this;
}

template<typename Tindices>
void SparseIJVMatrix<Tindices>::reserve(size_t size)
{
	_rowIndices.reserve(size);
	_columnIndices.reserve(size);
	_values.reserve(size);
}

template<typename Tindices>
void SparseIJVMatrix<Tindices>::transpose()
{
	_rowIndices.swap(_columnIndices);
	size_t tmp = _rows;
	_rows = _columns;
	_columns = tmp;
	sort();
}

template<typename Tindices>
void SparseIJVMatrix<Tindices>::sort(size_t begin, size_t end)
{
	if (begin >= end) {
		return;
	}

	size_t l = begin, h = end;

	do {
		while ((l < h) && compare(l, end) <= 0) { l++; } // lesser or equal
		while ((l < h) && compare(h, end) >= 0) { h--; } // higher of equal
		if (l < h) { swap(l, h); }
	} while (l < h);
	// swap pivot to the middle of the array
	swap(h, end);

	/* Sort smaller array first for less stack usage */
	if (l - begin < end - l) {
		sort(begin, l ? l - 1 : 0);
		sort(h + 1, end);
	} else {
		sort(h + 1, end);
		sort(begin, l ? l - 1 : 0);
	}
}

}

