#ifndef SPARSEIJVMATRIX_H_
#define SPARSEIJVMATRIX_H_

#include "matrix.h"
#include "denseMatrix.h"
#include "sparseDOKMatrix.h"
#include "sparseCSRMatrix.h"
#include "sparseVVPMatrix.h"

class DenseMatrix;
class SparseDOKMatrix;
class SparseCSRMatrix;
class SparseVVPMatrix;

class SparseIJVMatrix: public Matrix
{

public:

	SparseIJVMatrix(MatrixType type, MKL_INT rowsAndCols): Matrix(type, rowsAndCols, rowsAndCols) { };
	SparseIJVMatrix(MatrixType type, MKL_INT rows, MKL_INT cols): Matrix(type, rows, cols) { };
	SparseIJVMatrix(MKL_INT rows, MKL_INT cols): Matrix(Matrix::GENERAL, rows, cols) { };

	SparseIJVMatrix(const DenseMatrix &other);
	SparseIJVMatrix(const SparseDOKMatrix &other);
	SparseIJVMatrix(const SparseCSRMatrix &other);
	SparseIJVMatrix(SparseVVPMatrix &other);

	SparseIJVMatrix& operator=(const DenseMatrix &other)
	{
		SparseIJVMatrix tmp(other);
		assign(*this, tmp);
		return *this;
	}

	SparseIJVMatrix& operator=(const SparseDOKMatrix &other)
	{
		SparseIJVMatrix tmp(other);
		assign(*this, tmp);
		return *this;
	}

	SparseIJVMatrix& operator=(const SparseCSRMatrix &other)
	{
		SparseIJVMatrix tmp(other);
		assign(*this, tmp);
		return *this;
	}

	SparseIJVMatrix& operator=(SparseVVPMatrix &other)
	{
		SparseIJVMatrix tmp(other);
		assign(*this, tmp);
		return *this;
	}

	MKL_INT nonZeroValues() const
	{
		return _values.size();
	}

	double operator()(MKL_INT row, MKL_INT column) const
	{
		arrange(row, column);
		for(int i = 0; i < _rowIndices.size(); i++)
		{
			if (_rowIndices[i] == row && _columnIndices[i] == column) {
				return _values[i];
			}
		}
		return 0;
	}

	const MKL_INT* rowIndices() const
	{
		return &_rowIndices[0];
	}

	const MKL_INT* columnIndices() const
	{
		return &_columnIndices[0];
	}

	const double* values() const
	{
		return &_values[0];
	}

	MKL_INT* rowIndices()
	{
		return &_rowIndices[0];
	}

	MKL_INT* columnIndices()
	{
		return &_columnIndices[0];
	}

	double* values()
	{
		return &_values[0];
	}

protected:

	void makeTransposition();

private:

	static void assign(SparseIJVMatrix &m1, SparseIJVMatrix &m2)
	{
		Matrix::assign(m1, m2);
		m1._rowIndices.swap(m2._rowIndices);
		m1._columnIndices.swap(m2._columnIndices);
		m1._values.swap(m2._values);

	}

	// Sparse COO data
	std::vector<MKL_INT> _rowIndices;
	std::vector<MKL_INT> _columnIndices;
	std::vector<double> _values;
};

#endif /* SPARSEIJVMATRIX_H_ */
