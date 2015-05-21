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

#define IJVMatrixIndexing Matrix::ZeroBased

class SparseIJVMatrix: public Matrix
{

public:

	SparseIJVMatrix(): Matrix(IJVMatrixIndexing) { };
	SparseIJVMatrix(MKL_INT rows, MKL_INT columns): Matrix(rows, columns, IJVMatrixIndexing) { };

	SparseIJVMatrix(const DenseMatrix &other);
	SparseIJVMatrix(const SparseDOKMatrix &other);
	SparseIJVMatrix(const SparseCSRMatrix &other);
	SparseIJVMatrix(SparseVVPMatrix &other);

	SparseIJVMatrix& operator=(const DenseMatrix &other);
	SparseIJVMatrix& operator=(const SparseDOKMatrix &other);
	SparseIJVMatrix& operator=(const SparseCSRMatrix &other);
	SparseIJVMatrix& operator=(SparseVVPMatrix &other);

	void reserve(size_t size);
	void transpose();

	size_t nonZeroValues() const
	{
		return _values.size();
	}

	double operator()(size_t row, size_t column) const
	{
		return get(row, column);
	}

	double get(size_t row, size_t column) const
	{
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

	MKL_INT* rowIndices()
	{
		return &_rowIndices[0];
	}

	const MKL_INT* columnIndices() const
	{
		return &_columnIndices[0];
	}

	MKL_INT* columnIndices()
	{
		return &_columnIndices[0];
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

	double& operator()(size_t row, size_t column)
	{
		for(int i = 0; i < _rowIndices.size(); i++)
		{
			if (_rowIndices[i] == row && _columnIndices[i] == column) {
				return _values[i];
			}
		}
		_rowIndices.push_back(row);
		_columnIndices.push_back(column);
		_values.push_back(0);
		return _values.back();
	}

	void set(size_t row, size_t column, double value)
	{
		for(int i = 0; i < _rowIndices.size(); i++)
		{
			if (_rowIndices[i] == row && _columnIndices[i] == column) {
				_values[i] = value;
			}
		}
		_rowIndices.push_back(row);
		_columnIndices.push_back(column);
		_values.push_back(value);
	}

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
