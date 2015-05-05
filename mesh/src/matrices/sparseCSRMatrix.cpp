#include "sparseCSRMatrix.h"

SparseCSRMatrix::SparseCSRMatrix(const DenseMatrix &other): Matrix(other.type(), other.rows(), other.columns())
{
	MKL_INT nnz = other.nonZeroValues();
	_rowPtrs.resize(other.rows() + 1);
	_columnIndices.resize(nnz);
	_values.resize(nnz);

	int info;
	MKL_INT job[6] = {
		0,										// convert from dense to CSR
		0,										// zero based indexing of input matrix
		0,										// zero based indexing of output matrix
		(_type == Matrix::SYMETRIC)? 1: 2,		// full matrix of upper triangular
		nnz,									// number of non-zero values
		1										// generate full output
	};

	mkl_ddnscsr (
		job, &_rows, &_cols,
		const_cast<double*>(other.values()), &_cols,
		values(), columnIndices(), rowPtrs(),
		&info);
}

SparseCSRMatrix::SparseCSRMatrix(const SparseDOKMatrix &other): Matrix(other.type(), other.rows(), other.columns())
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

		std::fill(_rowPtrs.begin() + last_index, _rowPtrs.begin() + row->first + 1, nnz);
		last_index = row->first + 1;

		ColumnMap::const_iterator column;
		for(column = columns.begin(); column != columns.end(); ++column) {
			if (column->second != 0) {
				_columnIndices.push_back(column->first);
				_values.push_back(column->second);
				nnz++;
			}
		}
	}
	std::fill(_rowPtrs.begin() + last_index, _rowPtrs.end(), nnz);
}

SparseCSRMatrix::SparseCSRMatrix(const SparseIJVMatrix &other): Matrix(other.type(), other.rows(), other.columns())
{
	MKL_INT nnz = other.nonZeroValues();
	_rowPtrs.resize(other.rows() + 1);
	_columnIndices.resize(nnz);
	_values.resize(nnz);

	MKL_INT job[6] = {
		2, 			// IJV to sorted CSR
		0,			// zero based indexing
		0,			// zero based indexing
		0,			// without any meaning
		nnz,		// non-zero values
		0,			// fill all output arrays
	};

	MKL_INT info;

	mkl_dcsrcoo(
		job, &_rows,
		values(), columnIndices(), rowPtrs(), &nnz,
		const_cast<double*>(other.values()), const_cast<MKL_INT*>(other.rowIndices()), const_cast<MKL_INT*>(other.columnIndices()),
		&info);
}

SparseCSRMatrix::SparseCSRMatrix(SparseVVPMatrix &other): Matrix(other.type(), other.rows(), other.columns())
{
	other.shrink();
	MKL_INT nnz = other.nonZeroValues();
	_rowPtrs.reserve(other.rows() + 1);
	_columnIndices.reserve(nnz);
	_values.reserve(nnz);

	const VVP &values = other.values();

	_rowPtrs.push_back(0);
	for (size_t row = 0; row < _rows; row++) {
		for (size_t column = 0; column < values[row].size(); column++) {
			_columnIndices.push_back(values[row][column].first);
			_values.push_back(values[row][column].second);
		}
		_rowPtrs.push_back(_values.size());
	}
}

void SparseCSRMatrix::makeTransposition()
{
	std::vector<MKL_INT> colPtrs, rowIndices;
	std::vector<double> values;

	colPtrs.resize(_cols + 1);
	rowIndices.resize(_columnIndices.size());
	values.resize(_values.size());

	std::vector<int> colCounters(_cols, 0);
	for (size_t i = 0; i < _columnIndices.size(); i++) {
		colCounters[_columnIndices[i]]++;
	}

	colPtrs[0] = 0;
	for (size_t i = 0; i < _cols; i++) {
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

