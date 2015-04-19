#include "sparseDOKMatrix.h"

SparseDOKMatrix::SparseDOKMatrix(const DenseMatrix &other): EditableMatrix(other.type(), other.rows(), other.columns())
{
	for (size_t r = 0; r < other.rows(); r++) {
		for (size_t c = (_type == Matrix::SYMETRIC)? r: 0; c < other.columns(); c++) {
			if (other(r, c) != 0) {
				_values[r][c] = other(r, c);
			}
		}
	}
}


SparseDOKMatrix::SparseDOKMatrix(const SparseCSRMatrix &other): EditableMatrix(other.type(), other.rows(), other.columns())
{
	const MKL_INT *rowPrts = other.rowPtrs();
	const MKL_INT *columnIndices = other.columnIndices();
	const double *values = other.values();

	for (size_t r = 0; r < other.rows(); r++) {
		for (size_t c = rowPrts[r]; c < rowPrts[r + 1]; c++) {
			if (values[c] != 0) {
				_values[r][columnIndices[c]] = values[c];
			}
		}
	}
}

SparseDOKMatrix::SparseDOKMatrix(const SparseIJVMatrix &other): EditableMatrix(other.type(), other.rows(), other.columns())
{
	const MKL_INT *rowIndices = other.rowIndices();
	const MKL_INT *columnIndices = other.columnIndices();
	const double *values = other.values();

	for (size_t i = 0; i < other.nonZeroValues(); i++) {
		if (values[i] != 0) {
			_values[rowIndices[i]][columnIndices[i]] = values[i];
		}
	}
}

MKL_INT SparseDOKMatrix::removeZeros()
{
	MKL_INT nnz = 0;
	MatrixMap::iterator row;
	for (row = _values.begin(); row != _values.end(); ++row) {
		ColumnMap &columns = row->second;

		ColumnMap::iterator column = columns.begin();
		ColumnMap::iterator to_remove;
		while (column != columns.end()) {
			if (column->second == 0) {
				to_remove = column;
				++column;
				columns.erase(to_remove);
			} else {
				++column;
				++nnz;
			}
		}
	}
	return nnz;
}

MKL_INT SparseDOKMatrix::nonZeroValues() const
{
	MKL_INT nnz = 0;
	MatrixMap::const_iterator row;
	for (row = _values.begin(); row != _values.end(); ++row) {
		nnz += row->second.size();
	}
	return nnz;
}

void SparseDOKMatrix::makeTransposition()
{
	MatrixMap values;

	MatrixMap::const_iterator row;
	ColumnMap::const_iterator column;
	for (row = _values.begin(); row != _values.end(); ++row) {
		for (column = row->second.begin(); column != row->second.end(); ++column) {
			values[column->first][row->first] = column->second;
		}
	}

	_values.swap(values);
}



