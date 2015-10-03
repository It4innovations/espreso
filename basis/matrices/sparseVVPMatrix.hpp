#include "sparseVVPMatrix.h"

template<typename Tindices>
void SparseVVPMatrix<Tindices>::resize(size_t rows, size_t columns)
{
	_rows = rows;
	_columns = columns;
	_values.resize(rows);
}

template<typename Tindices>
void SparseVVPMatrix<Tindices>::shrink()
{
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
}

template<typename Tindices>
size_t SparseVVPMatrix<Tindices>::nonZeroValues() const
{
	size_t size = 0;
	for (size_t row = 0; row < _rows; row++) {
		size += _values[row].size();
	}
	return size;
}

template<typename Tindices>
void SparseVVPMatrix<Tindices>::transpose()
{
	//TODO: implement
}



