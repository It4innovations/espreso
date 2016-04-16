#include "sparseVVPMatrix.h"

namespace espreso {

template<typename Tindices>
std::ostream& operator<<(std::ostream& os, const SparseVVPMatrix<Tindices> &m)
{
	os << m.rows() << " " << m.columns() << " " << m.nonZeroValues() << "\n";

	for (size_t r = 0; r < m._values.size(); r++) {
		os << r << ": ";
		for (size_t c = 0; c < m._values[r].size(); c++) {
			os << "(" << m._values[r][c].first << ":" << m._values[r][c].second << ")";
		}
		os << "\n";
	}
	return os;
}

template<typename Tindices>
void SparseVVPMatrix<Tindices>::resize(size_t rows, size_t columns)
{
	_rows = rows;
	_columns = columns;
	_values.resize(rows);
}

static bool compare(const std::pair<eslocal, double> &a, const std::pair<eslocal, double> &b)
{
	return a.first < b.first;
}

template<typename Tindices>
void SparseVVPMatrix<Tindices>::shrink()
{
	size_t unique;
	for (size_t row = 0; row < _rows; row++) {
		std::sort(_values[row].begin(), _values[row].end(), compare);

		unique = 0;
		for (size_t j = 1; j < _values[row].size(); j++) {
			if (_values[row][j - 1].first == _values[row][j].first) {
				_values[row][unique].second += _values[row][j].second;
			} else {
				unique++;
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

}


