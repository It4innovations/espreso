#include "sparseCSRMatrix.h"

SparseCSRMatrix::SparseCSRMatrix(): Matrix(CSRMatrixIndexing)
{
	_rowPtrs.assign(2, _indexing);
}

SparseCSRMatrix::SparseCSRMatrix(size_t rows, size_t columns): Matrix(rows, columns, CSRMatrixIndexing)
{
	_rowPtrs.assign(rows + 1, _indexing);
}

SparseCSRMatrix::SparseCSRMatrix(const DenseMatrix &other): Matrix(other.rows(), other.columns(), CSRMatrixIndexing)
{
	MKL_INT nnz = other.nonZeroValues();
	MKL_INT rows = _rows;
	MKL_INT columns = _columns;
	_rowPtrs.resize(other.rows() + 1);
	_columnIndices.resize(nnz);
	_values.resize(nnz);

	int info;
	MKL_INT job[6] = {
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

SparseCSRMatrix::SparseCSRMatrix(const SparseDOKMatrix &other): Matrix(other.rows(), other.columns(), CSRMatrixIndexing)
{
	MKL_INT nnz = other.nonZeroValues();
	_rowPtrs.resize(other.rows() + 1);
	_columnIndices.reserve(nnz);
	_values.reserve(nnz);

	MKL_INT last_index = 0;

	const MatrixMap &dokValues = other.values();
	MatrixMap::const_iterator row;

	nnz = 0;
	for(row = dokValues.begin(); row != dokValues.end(); ++row) {
		const ColumnMap &columns = row->second;

		std::fill(_rowPtrs.begin() + last_index, _rowPtrs.begin() + row->first - other.indexing() + 1, nnz + _indexing);
		last_index = row->first - other.indexing() + 1;

		ColumnMap::const_iterator column;
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

SparseCSRMatrix::SparseCSRMatrix(const SparseIJVMatrix &other): Matrix(other.rows(), other.columns(), CSRMatrixIndexing)
{
	MKL_INT nnz = other.nonZeroValues();
	MKL_INT rows = _rows;
	_rowPtrs.resize(other.rows() + 1);
	_columnIndices.resize(nnz);
	_values.resize(nnz);

	MKL_INT job[6] = {
		2, 					// IJV to sorted CSR
		indexing(),			// indexing of CSR matrix
		other.indexing(),	// indexing of IJV matrix
		0,					// without any meaning
		nnz,				// non-zero values
		0,					// fill all output arrays
	};

	MKL_INT info;

	mkl_dcsrcoo(
		job, &rows,
		values(), columnIndices(), rowPtrs(), &nnz,
		const_cast<double*>(other.values()), const_cast<MKL_INT*>(other.rowIndices()), const_cast<MKL_INT*>(other.columnIndices()),
		&info);
}

SparseCSRMatrix::SparseCSRMatrix(SparseVVPMatrix &other): Matrix(other.rows(), other.columns(), CSRMatrixIndexing)
{
	other.shrink();
	MKL_INT nnz = other.nonZeroValues();
	_rowPtrs.reserve(other.rows() + 1);
	_columnIndices.reserve(nnz);
	_values.reserve(nnz);

	const VVP &values = other.values();

	_rowPtrs.push_back(_indexing);
	for (size_t row = 0; row < _rows; row++) {
		for (size_t column = 0; column < values[row].size(); column++) {
			_columnIndices.push_back(values[row][column].first + _indexing);
			_values.push_back(values[row][column].second);
		}
		_rowPtrs.push_back(_values.size() + _indexing);
	}
}

SparseCSRMatrix& SparseCSRMatrix::operator=(const DenseMatrix &other)
{
	SparseCSRMatrix tmp(other);
	assign(*this, tmp);
	return *this;
}

SparseCSRMatrix& SparseCSRMatrix::operator=(const SparseDOKMatrix &other)
{
	SparseCSRMatrix tmp(other);
	assign(*this, tmp);
	return *this;
}

SparseCSRMatrix& SparseCSRMatrix::operator=(const SparseIJVMatrix &other)
{
	SparseCSRMatrix tmp(other);
	assign(*this, tmp);
	return *this;
}

SparseCSRMatrix& SparseCSRMatrix::operator=(SparseVVPMatrix &other)
{
	SparseCSRMatrix tmp(other);
	assign(*this, tmp);
	return *this;
}

void SparseCSRMatrix::multiply(SparseCSRMatrix &A, SparseCSRMatrix &B, bool transposeA)
{
	if (_indexing == Matrix::ZeroBased) {
		std::cerr << "Multiplication of two CSR matrices with zero based indexing is not supported\n";
		exit(EXIT_FAILURE);
	}

	_rows = A.rows();
	_columns = B.columns();

	MKL_INT request = 0;
	MKL_INT sort = 8;	// C is sorted
	MKL_INT m = A.rows();
	MKL_INT n = A.columns();
	MKL_INT k = B.columns();
	MKL_INT nnz;		// not used
	MKL_INT info;

	_rowPtrs.resize(A.rows() + 1);

	request = 1;	// compute only rowPrts
	mkl_dcsrmultcsr(
		"n",
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

void SparseCSRMatrix::resize(size_t rows, size_t columns)
{
	_rows = rows;
	_columns = columns;
	_rowPtrs.resize(rows + 1, _rowPtrs.back());
}

void SparseCSRMatrix::transpose()
{
	std::vector<MKL_INT> colPtrs, rowIndices;
	std::vector<double> values;

	colPtrs.resize(_columns + 1);
	rowIndices.resize(_columnIndices.size());
	values.resize(_values.size());

	std::vector<MKL_INT> colCounters(_columns, 0);
	for (size_t i = 0; i < _columnIndices.size(); i++) {
		colCounters[_columnIndices[i]]++;
	}

	colPtrs[0] = 0;
	for (size_t i = 0; i < _columns; i++) {
		colPtrs[i + 1] = colPtrs[i] + colCounters[i];
	}
	std::fill(colCounters.begin(), colCounters.end(), 0);

	for (size_t r = 0; r < _rows; r++) {
		for (MKL_INT c = _rowPtrs[r]; c < _rowPtrs[r + 1]; c++) {
			MKL_INT index = colPtrs[_columnIndices[c]] + colCounters[_columnIndices[c]]++;
			rowIndices[index] = r;
			values[index] = _values[c];
		}
	}

	_rowPtrs.swap(colPtrs);
	_columnIndices.swap(rowIndices);
	_values.swap(values);
}

