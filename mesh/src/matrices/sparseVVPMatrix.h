#ifndef SPARSEVVPMATRIX_H_
#define SPARSEVVPMATRIX_H_

#include <algorithm>

#include "matrix.h"

typedef std::vector<std::vector<std::pair<MKL_INT, double> > > VVP;

class SparseVVPMatrix: public EditableMatrix
{

public:
	SparseVVPMatrix(MatrixType type, MKL_INT rowsAndCols)
		: EditableMatrix(type, rowsAndCols, rowsAndCols), _values(rowsAndCols), _shrunk(true) { };
	SparseVVPMatrix(MatrixType type, MKL_INT rows, MKL_INT cols)
		: EditableMatrix(type, rows, cols), _values(rows), _shrunk(true) { };
	SparseVVPMatrix(MKL_INT rows, MKL_INT cols)
		: EditableMatrix(Matrix::GENERAL, rows, cols), _values(rows), _shrunk(true) { };
	SparseVVPMatrix()
		: EditableMatrix(Matrix::GENERAL, 0, 0), _shrunk(true) { };

	double& operator()(MKL_INT row, MKL_INT column)
	{
		arrange(row, column);
		_shrunk = false;
		_values[row].push_back(std::pair<MKL_INT, double>(column, 0));
		return _values[row].back().second;
	}

	bool isShrunk()
	{
		return _shrunk;
	}

	void resize(MKL_INT rows, MKL_INT columns);
	void shrink();

	const VVP& values() const
	{
		return _values;
	}

	MKL_INT nonZeroValues() const;

	void makeTransposition();

private:

	double operator()(MKL_INT row, MKL_INT column) const
	{
		arrange(row, column);
		double value;
		for(int i = 0; i < _values[row].size(); i++)
		{
			if (_values[row][i].first == column) {
				value += _values[row][i].second;
			}
		}
		return value;
	}

	VVP _values;

	bool _shrunk;
};



#endif /* SPARSEVVPMATRIX_H_ */
