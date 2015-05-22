#ifndef DENSEMATRIX_H_
#define DENSEMATRIX_H_

#include <cstdio>
#include <algorithm>

#include "matrix.h"
#include "sparseDOKMatrix.h"
#include "sparseCSRMatrix.h"
#include "sparseIJVMatrix.h"

#define DenseMatrixIndexing Matrix::ZeroBased

class DenseMatrix: public Matrix
{

public:

	DenseMatrix(): Matrix(DenseMatrixIndexing) {};
	DenseMatrix(size_t rows, size_t columns)
		: Matrix(rows, columns, DenseMatrixIndexing), _values(rows * columns, 0) {};

	DenseMatrix(const SparseDOKMatrix &other);
	DenseMatrix(const SparseCSRMatrix &other);
	DenseMatrix(const SparseIJVMatrix &other);

	DenseMatrix& operator=(const SparseDOKMatrix &other);
	DenseMatrix& operator=(const SparseCSRMatrix &other);
	DenseMatrix& operator=(const SparseIJVMatrix &other);

	void multiply(DenseMatrix &A, DenseMatrix &B, double alfa = 1, double beta = 0, bool transposeA = false, bool transposeB = false);

	void resize(size_t rows, size_t columns);
	void transpose();

	size_t nonZeroValues() const
	{
		return std::count_if(_values.begin(), _values.end(), NonZeroValue());
	}

	double operator()(size_t row, size_t column) const
	{
		return _values[row * _columns + column];
	};

	double& operator()(size_t row, size_t column)
	{
		return _values[row * _columns + column];
	}

	double get(size_t row, size_t column) const
	{
		return _values[row * _columns + column];
	}

	void set(size_t row, size_t column, double value)
	{
		_values[row * _columns + column] = value;
	}

	const double* values() const
	{
		return &_values[0];
	}

	double* values()
	{
		return &_values[0];
	}

private:

	static void assign(DenseMatrix &m1, DenseMatrix &m2)
	{
		Matrix::assign(m1, m2);
		m1._values.swap(m2._values);
	}

	std::vector<double> _values;
};


#endif /* DENSEMATRIX_H_ */
