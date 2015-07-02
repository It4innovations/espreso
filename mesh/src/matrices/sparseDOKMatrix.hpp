#include "sparseDOKMatrix.h"

template<typename Tindices>
size_t SparseDOKMatrix<Tindices>::nonZeroValues() const
{
	size_t nnz = 0;
	typename std::map<Tindices, std::map<Tindices, double> >::const_iterator row;
	for (row = _values.begin(); row != _values.end(); ++row) {
		nnz += row->second.size();
	}
	return nnz;
}

template<typename Tindices>
void SparseDOKMatrix<Tindices>::transpose()
{
	typename std::map<Tindices, std::map<Tindices, double> > values;

	typename std::map<Tindices, std::map<Tindices, double> >::const_iterator row;
	typename std::map<Tindices, double>::const_iterator column;
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



