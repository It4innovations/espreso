#ifndef SPARSEVVPMATRIX_H_
#define SPARSEVVPMATRIX_H_

#include <algorithm>

#include "matrix.h"

#define VVPMatrixIndexing Matrix::ZeroBased

template<typename Tindices>
using VVP = std::vector<std::vector<std::pair<Tindices, double> > >;

template<typename Tindices>
class SparseVVPMatrix: public Matrix
{

public:

	SparseVVPMatrix(): Matrix(VVPMatrixIndexing) {};
	SparseVVPMatrix(size_t rows, size_t columns): Matrix(rows, columns, VVPMatrixIndexing), _values(rows) {};

	void shrink();
	void resize(size_t rows, size_t columns);
	void transpose();
	size_t nonZeroValues() const;

	double& operator()(size_t row, size_t column)
	{
		if(_values.size() <= row) {
			_values.resize(row + 1);
			_rows = row + 1;
		}
		if (column >= _columns) {
			_columns = column + 1;
		}
		_values[row].push_back(std::pair<Tindices, double>(column, 0));
		return _values[row].back().second;
	}

	void set(size_t row, size_t column, double value)
	{
		if (Matrix::nonZero(value)) {
			if(_values.size() <= row) {
				_values.resize(row + 1);
				_rows = row + 1;
			}
			if (column >= _columns) {
				_columns = column + 1;
			}
			_values[row].push_back(std::pair<Tindices, double>(column, value));
		}
	}

	const VVP<Tindices>& values() const
	{
		return _values;
	}


private:

	double operator()(size_t row, size_t column) const
	{
		return get(row, column);
	}

	double get(size_t row, size_t column) const
	{
		double value = 0;
		for (size_t i = 0; i < _values[row].size(); i++) {
			if (_values[row][i].first == column) {
				value += _values[row][i].second;
			}
		}
		return value;
	}

	VVP<Tindices> _values;
};

#include "sparseVVPMatrix.hpp"



#endif /* SPARSEVVPMATRIX_H_ */
