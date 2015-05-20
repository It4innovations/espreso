#ifndef SPARSEVVPMATRIX_H_
#define SPARSEVVPMATRIX_H_

#include <algorithm>

#include "matrix.h"

typedef std::vector<std::vector<std::pair<size_t, double> > > VVP;

class SparseVVPMatrix: public Matrix
{

public:

	SparseVVPMatrix() {};
	SparseVVPMatrix(size_t rows, size_t columns): Matrix(rows, columns), _values(rows) {};

	void shrink();
	void resize(size_t rows, size_t columns);
	void transpose();
	size_t nonZeroValues() const;

	double& operator()(size_t row, size_t column)
	{
		_values[row].push_back(std::pair<size_t, double>(column, 0));
		return _values[row].back().second;
	}

	void set(size_t row, size_t column, double value)
	{
		_values[row].push_back(std::pair<size_t, double>(column, value));
	}

	const VVP& values() const
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

	VVP _values;
};



#endif /* SPARSEVVPMATRIX_H_ */
