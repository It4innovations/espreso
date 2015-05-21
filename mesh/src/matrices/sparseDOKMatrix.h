#ifndef SPARSEDOKMATRIX_H_
#define SPARSEDOKMATRIX_H_

#include <map>

#include "matrix.h"
#include "denseMatrix.h"
#include "sparseCSRMatrix.h"
#include "sparseIJVMatrix.h"

#define DOKMatrixIndexing Matrix::ZeroBased

typedef std::map<size_t, std::map<size_t, double> > MatrixMap;
typedef std::map<size_t, double> ColumnMap;

class SparseDOKMatrix: public Matrix
{

public:

	SparseDOKMatrix(): Matrix(DOKMatrixIndexing) {};
	SparseDOKMatrix(size_t rows, size_t columns): Matrix(rows, columns, DOKMatrixIndexing) {};

	void transpose();
	size_t nonZeroValues() const;


	double operator()(size_t row, size_t column) const
	{
		return get(row, column);
	};

	double& operator () (size_t row, size_t column)
	{
		return _values[row][column];
	}

	double get(size_t row, size_t column) const
	{
		MatrixMap::const_iterator row_it = _values.find(row);
		if (row_it == _values.end()) {
			return 0;
		}
		ColumnMap::const_iterator column_it = row_it->second.find(column);
		if (column_it == row_it->second.end()) {
			return 0;
		}
		return column_it->second;
	}
	void set(size_t row, size_t column, double value)
	{
		if (Matrix::nonZero(value)) {
			_values[row][column] = value;
		}
	}

	const MatrixMap& values() const
	{
		return _values;
	}

	MatrixMap& values()
	{
		return _values;
	}

private:

	static void assign(SparseDOKMatrix &m1, SparseDOKMatrix &m2)
	{
		Matrix::assign(m1, m2);
		m1._values.swap(m2._values);
	}

	MatrixMap _values;
};

#endif /* SPARSEDOKVMATRIX_H_ */
