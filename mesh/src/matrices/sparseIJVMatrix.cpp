#include "sparseIJVMatrix.h"

std::ostream& operator<<(std::ostream& os, const SparseIJVMatrix &m)
{
	os << m.rows() << " " << m.columns() << " " << m.nonZeroValues() << "\n";

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
	//TODO: po swap by mela byt matice znovu setridena
	_rowIndices.swap(_columnIndices);
	size_t tmp = _rows;
	_rows = _columns;
	_columns = tmp;
	sort();
}

#define compare(a, b, i, j) (a[i] != a[j] ? a[i] - a[j] : b[i] - b[j])

void SparseIJVMatrix::sort(size_t begin, size_t end)
{
	size_t l, h, p;
	eslocal p1, p2, t;
	double p3, td;

	if (begin >= end) {
		return;
	}

	l = begin;
	h = end;
	p = end;

	p1 = _columnIndices[p];
	p2 = _rowIndices[p];
	p3 = _values[p];

	do {
		while ((l < h) && compare(_rowIndices, _columnIndices, l, p) <= 0) { l++; }
		while ((h > l) && compare(_rowIndices, _columnIndices, h, p) >= 0) { h--; }
		if (l < h) {
			t = _columnIndices[l];
			_columnIndices[l] = _columnIndices[h];
			_columnIndices[h] = t;

			t = _rowIndices[l];
			_rowIndices[l] = _rowIndices[h];
			_rowIndices[h] = t;

			td = _values[l];
			_values[l] = _values[h];
			_values[h] = td;
		}
	} while (l < h);

	_columnIndices[p] = _columnIndices[l];
	_columnIndices[l] = p1;

	_rowIndices[p] = _rowIndices[l];
	_rowIndices[l] = p2;

	_values[p] = _values[l];
	_values[l] = p3;

	/* Sort smaller array first for less stack usage */
	if (l - begin < end - l) {
		sort(begin, l - 1);
		sort(l + 1, end);
	} else {
		sort(l + 1, end);
		sort(begin, l - 1);
	}
}
