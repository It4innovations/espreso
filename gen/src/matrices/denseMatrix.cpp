#include "denseMatrix.h"

DenseMatrix::DenseMatrix(MatrixType type, MKL_INT rowsAndCols): EditableMatrix(type, rowsAndCols, rowsAndCols)
{
	std::vector<double>(rows() * columns(), 0).swap(_values);
}

DenseMatrix::DenseMatrix(const SparseDOKMatrix &other): EditableMatrix(other.type(), other.rows(), other.columns())
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

DenseMatrix::DenseMatrix(const SparseCSRMatrix &other): EditableMatrix(other.type(), other.rows(), other.columns())
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

DenseMatrix::DenseMatrix(const SparseIJVMatrix &other): EditableMatrix(other.type(), other.rows(), other.columns())
{
	std::vector<double>(rows() * columns(), 0).swap(_values);

	const MKL_INT *rowIndices = other.rowIndices();
	const MKL_INT *columnIndices = other.columnIndices();
	const double *values = other.values();

	for (MKL_INT i = 0; i < other.nonZeroValues(); i++) {
		this->operator ()(rowIndices[i], columnIndices[i]) = values[i];
	}
}

void DenseMatrix::makeTransposition()
{
	std::vector<double> values(_rows * _cols);
	for (MKL_INT r = 0; r < _rows; r++) {
		for (MKL_INT c = 0; c < _cols; c++){
			values[c * _rows + r] = _values[r * _cols + c];
		}
	}

	_values.swap(values);
}
