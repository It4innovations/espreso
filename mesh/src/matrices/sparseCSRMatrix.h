#ifndef SPARSECSRMATRIX_H_
#define SPARSECSRMATRIX_H_

#include "matrix.h"
#include "sparseDOKMatrix.h"
#include "sparseIJVMatrix.h"
#include "sparseVVPMatrix.h"

class DenseMatrix;
class SparseDOKMatrix;
class SparseIJVMatrix;
class SparseVVPMatrix;

class SparseCSRMatrix: public Matrix
{

public:

	SparseCSRMatrix() { };
	SparseCSRMatrix(size_t rows, size_t columns): Matrix(rows, columns) { };

	SparseCSRMatrix(const DenseMatrix &other);
	SparseCSRMatrix(const SparseDOKMatrix &other);
	SparseCSRMatrix(const SparseIJVMatrix &other);
	SparseCSRMatrix(SparseVVPMatrix &other);

	SparseCSRMatrix& operator=(const DenseMatrix &other);
	SparseCSRMatrix& operator=(const SparseDOKMatrix &other);
	SparseCSRMatrix& operator=(const SparseIJVMatrix &other);
	SparseCSRMatrix& operator=(SparseVVPMatrix &other);

	void resize(size_t rows, size_t values);
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
		for(int i = _rowPtrs[row]; i < _rowPtrs[row + 1]; i++) {
			if (_columnIndices[i] == column) {
				return _values[i];
			}
		}
		return 0;
	}

	const double* values() const
	{
		return &_values[0];
	}

	double* values()
	{
		return &_values[0];
	}

	const MKL_INT* rowPtrs() const
	{
		return &_rowPtrs[0];
	}

	MKL_INT* rowPtrs()
	{
		return &_rowPtrs[0];
	}

	const MKL_INT* columnIndices() const
	{
		return &_columnIndices[0];
	}

	MKL_INT* columnIndices()
	{
		return &_columnIndices[0];
	}

private:

	double& operator()(size_t row, size_t column)
	{
		std::cerr << "It is not possible to insert to CRS matrix.\n";
		exit(EXIT_FAILURE);
	}

	void set(size_t row, size_t column, double value)
	{
		std::cerr << "It is not possible to insert to CRS matrix.\n";
		exit(EXIT_FAILURE);
	}

	static void assign(SparseCSRMatrix &m1, SparseCSRMatrix &m2)
	{
		Matrix::assign(m1, m2);
		m1._rowPtrs.swap(m2._rowPtrs);
		m1._columnIndices.swap(m2._columnIndices);
		m1._values.swap(m2._values);
	}

	// Sparse CSR data
	std::vector<MKL_INT> _rowPtrs;
	std::vector<MKL_INT> _columnIndices;
	std::vector<double> _values;

};

#endif /* SPARSEIJVMATRIX_H_ */
