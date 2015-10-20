#ifndef SPARSEDOKMATRIX_H_
#define SPARSEDOKMATRIX_H_

#include <map>

#include "matrix.h"
#include "denseMatrix.h"
#include "sparseCSRMatrix.h"
#include "sparseIJVMatrix.h"

#define DOKMatrixIndexing Matrix::ZeroBased

template<typename Tindices>
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
	}

	double& operator () (size_t row, size_t column)
	{
		return _values[row + _indexing][column + _indexing];
	}

	double get(size_t row, size_t column) const
	{
		typename std::map<Tindices, std::map<Tindices, double> >::const_iterator row_it = _values.find(row + _indexing);
		if (row_it == _values.end()) {
			return 0;
		}
		typename std::map<Tindices, double>::const_iterator column_it = row_it->second.find(column + _indexing);
		if (column_it == row_it->second.end()) {
			return 0;
		}
		return column_it->second;
	}
	void set(size_t row, size_t column, double value)
	{
		if (Matrix::nonZero(value)) {
			_values[row + _indexing][column + _indexing] = value;
		}
	}

	const std::map<Tindices, std::map<Tindices, double> >& values() const
	{
		return _values;
	}

	std::map<Tindices, std::map<Tindices, double> >& values()
	{
		return _values;
	}

private:

	static void assign(SparseDOKMatrix<Tindices> &m1, SparseDOKMatrix<Tindices> &m2)
	{
		Matrix::assign(m1, m2);
		m1._values.swap(m2._values);
	}

	std::map<Tindices, std::map<Tindices, double> > _values;
};

#include "sparseDOKMatrix.hpp"

#endif /* SPARSEDOKVMATRIX_H_ */
