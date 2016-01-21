#include "sparseCSRMatrix.h"


template<typename Tindices>
std::ostream& operator<<(std::ostream& os, const SparseCSRMatrix<Tindices> &m)
{
	SparseIJVMatrix<Tindices> ijv = m;
	os << ijv;
//	os << m.rows() << " " << m.columns() << " " << m.nonZeroValues() << "\n";
//
//	os << m._rowPtrs << "\n";
//	os << m._columnIndices << "\n";
//	os << m._values << "\n";
	return os;
}

template<typename Tindices>
SparseCSRMatrix<Tindices>::SparseCSRMatrix(): Matrix(CSRMatrixIndexing)
{
	_rowPtrs.assign(2, _indexing);
}

template<typename Tindices>
SparseCSRMatrix<Tindices>::SparseCSRMatrix(size_t rows, size_t columns): Matrix(rows, columns, CSRMatrixIndexing)
{
	_rowPtrs.assign(rows + 1, _indexing);
}

template<typename Tindices>
SparseCSRMatrix<Tindices>::SparseCSRMatrix(const DenseMatrix &other): Matrix(other.rows(), other.columns(), CSRMatrixIndexing)
{
	Tindices nnz = other.nonZeroValues();
	Tindices rows = _rows;
	Tindices columns = _columns;
	_rowPtrs.resize(other.rows() + 1);
	_columnIndices.resize(nnz);
	_values.resize(nnz);

	Tindices info;
	Tindices job[6] = {
		0,					// convert from dense to CSR
		other.indexing(),	// indexing of dense matrix
		indexing(),			// indexing of CSR matrix
		2,					// full matrix
		nnz,				// number of non-zero values
		1					// generate full output
	};

	mkl_ddnscsr (
		job, &rows, &columns,
		const_cast<double*>(other.values()), &columns,
		values(), columnIndices(), rowPtrs(),
		&info);
}

template<typename Tindices>
SparseCSRMatrix<Tindices>::SparseCSRMatrix(const SparseDOKMatrix<Tindices> &other): Matrix(other.rows(), other.columns(), CSRMatrixIndexing)
{
	size_t nnz = other.nonZeroValues();
	_rowPtrs.resize(other.rows() + 1);
	_columnIndices.reserve(nnz);
	_values.reserve(nnz);

	Tindices last_index = 0;

	const std::map<Tindices, std::map<Tindices, double> > &dokValues = other.values();
	typename std::map<Tindices, std::map<Tindices, double> >::const_iterator row;

	nnz = 0;
	for(row = dokValues.begin(); row != dokValues.end(); ++row) {
		const typename std::map<Tindices, double> &columns = row->second;

		std::fill(_rowPtrs.begin() + last_index, _rowPtrs.begin() + row->first - other.indexing() + 1, nnz + _indexing);
		last_index = row->first - other.indexing() + 1;

		typename std::map<Tindices, double>::const_iterator column;
		for(column = columns.begin(); column != columns.end(); ++column) {
			if (column->second != 0) {
				_columnIndices.push_back(column->first - other.indexing() + _indexing);
				_values.push_back(column->second);
				nnz++;
			}
		}
	}
	std::fill(_rowPtrs.begin() + last_index, _rowPtrs.end(), nnz + _indexing);
}

template<typename Tindices>
SparseCSRMatrix<Tindices>::SparseCSRMatrix(const SparseIJVMatrix<Tindices> &other): Matrix(other.rows(), other.columns(), CSRMatrixIndexing)
{
	Tindices nnz = other.nonZeroValues();
	Tindices rows = _rows;
	_rowPtrs.resize(other.rows() + 1);
	_columnIndices.resize(nnz);
	_values.resize(nnz);

	Tindices job[6] = {
		2, 					// IJV to sorted CSR
		indexing(),			// indexing of CSR matrix
		other.indexing(),	// indexing of IJV matrix
		0,					// without any meaning
		nnz,				// non-zero values
		0,					// fill all output arrays
	};

	Tindices info;

	mkl_dcsrcoo(
		job, &rows,
		values(), columnIndices(), rowPtrs(), &nnz,
		const_cast<double*>(other.values()), const_cast<eslocal*>(other.rowIndices()), const_cast<eslocal*>(other.columnIndices()),
		&info);
}

template<typename Tindices>
SparseCSRMatrix<Tindices>::SparseCSRMatrix(SparseVVPMatrix<Tindices> &other): Matrix(other.rows(), other.columns(), CSRMatrixIndexing)
{
	other.shrink();
	size_t nnz = other.nonZeroValues();
	_rowPtrs.reserve(other.rows() + 1);
	_columnIndices.reserve(nnz);
	_values.reserve(nnz);

	const VVP<Tindices> &values = other.values();

	_rowPtrs.push_back(_indexing);
	for (size_t row = 0; row < _rows; row++) {
		for (size_t column = 0; column < values[row].size(); column++) {
			_columnIndices.push_back(values[row][column].first + _indexing);
			_values.push_back(values[row][column].second);
		}
		_rowPtrs.push_back(_values.size() + _indexing);
	}
}

template<typename Tindices>
SparseCSRMatrix<Tindices>& SparseCSRMatrix<Tindices>::operator=(const DenseMatrix &other)
{
	SparseCSRMatrix<Tindices> tmp(other);
	assign(*this, tmp);
	return *this;
}

template<typename Tindices>
SparseCSRMatrix<Tindices>& SparseCSRMatrix<Tindices>::operator=(const SparseDOKMatrix<Tindices> &other)
{
	SparseCSRMatrix<Tindices> tmp(other);
	assign(*this, tmp);
	return *this;
}

template<typename Tindices>
SparseCSRMatrix<Tindices>& SparseCSRMatrix<Tindices>::operator=(const SparseIJVMatrix<Tindices> &other)
{
	SparseCSRMatrix<Tindices> tmp(other);
	assign(*this, tmp);
	return *this;
}

template<typename Tindices>
SparseCSRMatrix<Tindices>& SparseCSRMatrix<Tindices>::operator=(SparseVVPMatrix<Tindices> &other)
{
	SparseCSRMatrix<Tindices> tmp(other);
	assign(*this, tmp);
	return *this;
}

template<typename Tindices>
void SparseCSRMatrix<Tindices>::multiply(SparseCSRMatrix<Tindices> &A, SparseCSRMatrix<Tindices> &B, bool transposeA)
{
	if (_indexing == Matrix::ZeroBased) {
		std::cerr << "Multiplication of two CSR matrices with zero based indexing is not supported\n";
		exit(EXIT_FAILURE);
	}
	if (transposeA) {
		if (A.rows() != B.rows()) {
			std::cerr << "Matrix multiplication: matrices have incorrect dimensions.\n";
			exit(EXIT_FAILURE);
		}
	} else {
		if (A.columns() != B.rows()) {
			std::cerr << "Matrix multiplication: matrices have incorrect dimensions.\n";
			exit(EXIT_FAILURE);
		}
	}

	Tindices request = 0;
	Tindices sort = 8;	// C is sorted
	Tindices m = A.rows();
	Tindices n = A.columns();
	Tindices k = B.columns();
	Tindices nnz;		// not used
	Tindices info;

	_rows = transposeA ? A.columns() : A.rows();;
	_columns = k;
	_rowPtrs.resize(_rows + 1);

	request = 1;	// compute only rowPrts
	mkl_dcsrmultcsr(
		transposeA ? "t" : "n",
		&request,
		&sort,
		&m, &n, &k,
		A.values(), A.columnIndices(), A.rowPtrs(),
		B.values(), B.columnIndices(), B.rowPtrs(),
		values(), columnIndices(), rowPtrs(),
		&nnz,
		&info);

	_columnIndices.resize(_rowPtrs.back() - 1);
	_values.resize(_rowPtrs.back() - 1);

	request = 2;	// compute the rest of the matrix
	mkl_dcsrmultcsr(
		transposeA ? "t" : "n",
		&request,
		&sort,
		&m, &n, &k,
		A.values(), A.columnIndices(), A.rowPtrs(),
		B.values(), B.columnIndices(), B.rowPtrs(),
		values(), columnIndices(), rowPtrs(),
		&nnz,
		&info);
}

template<typename Tindices>
void SparseCSRMatrix<Tindices>::resize(size_t rows, size_t columns)
{
	_rows = rows;
	_columns = columns;
	_rowPtrs.resize(rows + 1, _rowPtrs.back());
}

template<typename Tindices>
void SparseCSRMatrix<Tindices>::transpose()
{
	Tindices job[6] = {
			0,		// CSR to CSC
			_indexing,
			_indexing,
			0,
			0,
			1
	};

	size_t size;
	if (_rows < _columns) {
		size = _columns;
		_rowPtrs.resize(size + 1, _rowPtrs.back());
	} else {
		size = _rows;
	}

	std::vector<Tindices> colPtrs(size + 1);
	std::vector<Tindices> rowIndices(_columnIndices.size());
	std::vector<double> vals(_values.size());

	Tindices n = size;
	Tindices info;

	mkl_dcsrcsc(
			job, &n,
			values(), columnIndices(), rowPtrs(),
			&vals[0], &rowIndices[0], &colPtrs[0],
			&info);

	colPtrs.resize(_columns + 1);
	_rowPtrs.swap(colPtrs);
	_columnIndices.swap(rowIndices);
	_values.swap(vals);
	size_t tmp = _rows;
	_rows = _columns;
	_columns = tmp;
}

