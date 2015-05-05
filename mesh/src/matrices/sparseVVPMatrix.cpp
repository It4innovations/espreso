#include "sparseVVPMatrix.h"

void SparseVVPMatrix::resize(MKL_INT rows, MKL_INT columns)
{
	_rows = rows;
	_cols = columns;
	_values.resize(rows);
}

void SparseVVPMatrix::shrink()
{
	if (_shrunk) {
		return;
	}
	size_t unique;
	for (size_t row = 0; row < _rows; row++) {
		std::sort(_values[row].begin(), _values[row].end());

		unique = 0;
		for (size_t j = 1; j < _values[row].size(); j++) {
			if (_values[row][j - 1].first == _values[row][j].first) {
				_values[row][unique].second += _values[row][j].second;
			} else {
				if (_values[row][unique].second != 0) {
					unique++;
				}
				_values[row][unique] = _values[row][j];
			}
		}
		_values[row].resize(unique + 1);
	}
	_shrunk = true;
}

MKL_INT SparseVVPMatrix::nonZeroValues() const
{
	MKL_INT size = 0;
	for (size_t row = 0; row < _rows; row++) {
		size += _values[row].size();
	}
	return size;
}

void SparseVVPMatrix::makeTransposition()
{
	//TODO: implement
}



