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

#define CSRMatrixIndexing Matrix::OneBased

class SparseCSRMatrix: public Matrix
{

public:

	SparseCSRMatrix();
	SparseCSRMatrix(size_t rows, size_t columns);

	SparseCSRMatrix(const DenseMatrix &other);
	SparseCSRMatrix(const SparseDOKMatrix &other);
	SparseCSRMatrix(const SparseIJVMatrix &other);
	SparseCSRMatrix(SparseVVPMatrix &other);

	SparseCSRMatrix& operator=(const DenseMatrix &other);
	SparseCSRMatrix& operator=(const SparseDOKMatrix &other);
	SparseCSRMatrix& operator=(const SparseIJVMatrix &other);
	SparseCSRMatrix& operator=(SparseVVPMatrix &other);

	void multiply(SparseCSRMatrix &A, SparseCSRMatrix &B, bool transposeA = false);

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
		for(eslocal i = _rowPtrs[row]; i < _rowPtrs[row + 1]; i++) {
			if (_columnIndices[i - _indexing] == column + _indexing) {
				return _values[i - _indexing];
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

	const eslocal* rowPtrs() const
	{
		return &_rowPtrs[0];
	}

	eslocal* rowPtrs()
	{
		return &_rowPtrs[0];
	}

	const eslocal* columnIndices() const
	{
		return &_columnIndices[0];
	}

	eslocal* columnIndices()
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
	std::vector<eslocal> _rowPtrs;
	std::vector<eslocal> _columnIndices;
	std::vector<double> _values;

};

#endif /* SPARSEIJVMATRIX_H_ */
