#ifndef SPARSECSRMATRIX_H_
#define SPARSECSRMATRIX_H_

#include "matrix.h"
#include "sparseDOKMatrix.h"
#include "sparseIJVMatrix.h"
#include "sparseVVPMatrix.h"

namespace espreso {

class DenseMatrix;
template<typename Tindices> class SparseDOKMatrix;
template<typename Tindices> class SparseIJVMatrix;
template<typename Tindices> class SparseVVPMatrix;

#define CSRMatrixIndexing Matrix::OneBased

template<typename Tindices>
class SparseCSRMatrix: public Matrix
{

public:

	template<typename Tindices>
	friend std::ostream& operator<<(std::ostream& os, const SparseCSRMatrix<Tindices> &m);

	SparseCSRMatrix();
	SparseCSRMatrix(size_t rows, size_t columns);

	SparseCSRMatrix(const DenseMatrix &other);
	SparseCSRMatrix(const SparseDOKMatrix<Tindices> &other);
	SparseCSRMatrix(const SparseIJVMatrix<Tindices> &other);
	SparseCSRMatrix(SparseVVPMatrix<Tindices> &other);

	SparseCSRMatrix& operator=(const DenseMatrix &other);
	SparseCSRMatrix<Tindices>& operator=(const SparseDOKMatrix<Tindices> &other);
	SparseCSRMatrix<Tindices>& operator=(const SparseIJVMatrix<Tindices> &other);
	SparseCSRMatrix<Tindices>& operator=(SparseVVPMatrix<Tindices> &other);

	void multiply(SparseCSRMatrix<Tindices> &A, SparseCSRMatrix<Tindices> &B, bool transposeA = false);

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

	const Tindices* rowPtrs() const
	{
		return &_rowPtrs[0];
	}

	Tindices* rowPtrs()
	{
		return &_rowPtrs[0];
	}

	const Tindices* columnIndices() const
	{
		return &_columnIndices[0];
	}

	Tindices* columnIndices()
	{
		return &_columnIndices[0];
	}

	double norm() const
	{
		double n = 0;
		for (size_t i = 0; i < _values.size(); i++) {
				n += _values[i] * _values[i];
		}
		return sqrt(n);
	}

private:

	double& operator()(size_t row, size_t column)
	{
		ESINFO(ERROR) << "It is not possible to insert to CRS matrix.";
		exit(EXIT_FAILURE);
	}

	void set(size_t row, size_t column, double value)
	{
		ESINFO(ERROR) << "It is not possible to insert to CRS matrix.";
	}

	static void assign(SparseCSRMatrix<Tindices> &m1, SparseCSRMatrix<Tindices> &m2)
	{
		Matrix::assign(m1, m2);
		m1._rowPtrs.swap(m2._rowPtrs);
		m1._columnIndices.swap(m2._columnIndices);
		m1._values.swap(m2._values);
	}

	// Sparse CSR data
	std::vector<Tindices> _rowPtrs;
	std::vector<Tindices> _columnIndices;
	std::vector<double> _values;
};

}

#include "sparseCSRMatrix.hpp"

#endif /* SPARSEIJVMATRIX_H_ */
