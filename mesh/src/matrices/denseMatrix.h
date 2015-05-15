#ifndef DENSEMATRIX_H_
#define DENSEMATRIX_H_

#include <cstdio>
#include <algorithm>

#include "matrix.h"
#include "sparseDOKMatrix.h"
#include "sparseCSRMatrix.h"
#include "sparseIJVMatrix.h"

struct NonZero
{
	bool operator()(double value)
	{
		return value != 0;
	}
};

class DenseMatrix: public EditableMatrix
{

public:

	DenseMatrix(MatrixType type, MKL_INT rowsAndCols);
	DenseMatrix(MKL_INT rows, MKL_INT cols):
		EditableMatrix(Matrix::GENERAL, rows, cols), _values(rows * cols, 0) {};

	DenseMatrix(const SparseDOKMatrix &other);
	DenseMatrix(const SparseCSRMatrix &other);
	DenseMatrix(const SparseIJVMatrix &other);

	DenseMatrix& operator=(const SparseDOKMatrix &other)
	{
		DenseMatrix tmp(other);
		assign(*this, tmp);
		return *this;
	}

	DenseMatrix& operator=(const SparseCSRMatrix &other)
	{
		DenseMatrix tmp(other);
		assign(*this, tmp);
		return *this;
	}

	DenseMatrix& operator=(const SparseIJVMatrix &other)
	{
		DenseMatrix tmp(other);
		assign(*this, tmp);
		return *this;
	}

	void resize(MKL_INT rows, MKL_INT columns);

	MKL_INT nonZeroValues() const
	{
		return std::count_if(_values.begin(), _values.end(), NonZero());
	}

	double operator()(MKL_INT row, MKL_INT column) const
	{
		arrange(row, column);
		return _values[row * _cols + column];
	};

	double& operator()(MKL_INT row, MKL_INT column)
	{
		arrange(row, column);
		return _values[row * _cols + column];
	}

	const double* values() const
	{
		return &_values[0];
	}

	double* values()
	{
		return &_values[0];
	}

protected:

	void makeTransposition();

private:

	static void assign(DenseMatrix &m1, DenseMatrix &m2)
	{
		Matrix::assign(m1, m2);
		m1._values.swap(m2._values);
	}

	std::vector<double> _values;
};


#endif /* DENSEMATRIX_H_ */
