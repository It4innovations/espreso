#include "sparseDOKMatrix.h"

size_t SparseDOKMatrix::nonZeroValues() const
{
	size_t nnz = 0;
	MatrixMap::const_iterator row;
	for (row = _values.begin(); row != _values.end(); ++row) {
		nnz += row->second.size();
	}
	return nnz;
}

void SparseDOKMatrix::transpose()
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
	size_t tmp = _rows;
	_rows = _columns;
	_columns = tmp;
}



