#ifndef SPARSECSRMATRIX_H_
#define SPARSECSRMATRIX_H_

#include "matrix.h"
#include "sparseDOKMatrix.h"
#include "sparseIJVMatrix.h"

class DenseMatrix;
class SparseDOKMatrix;
class SparseIJVMatrix;

class SparseCSRMatrix: public Matrix
{

public:

	SparseCSRMatrix(MatrixType type, MKL_INT rowsAndCols): Matrix(type, rowsAndCols, rowsAndCols) { };
	SparseCSRMatrix(MatrixType type, MKL_INT rows, MKL_INT cols): Matrix(type, rows, cols) { };
	SparseCSRMatrix(MKL_INT rows, MKL_INT cols): Matrix(Matrix::GENERAL, rows, cols) { };

	SparseCSRMatrix(const DenseMatrix &other);
	SparseCSRMatrix(const SparseDOKMatrix &other);
	SparseCSRMatrix(const SparseIJVMatrix &other);

	SparseCSRMatrix& operator=(const DenseMatrix &other)
	{
		SparseCSRMatrix tmp(other);
		assign(*this, tmp);
		return *this;
	}

	SparseCSRMatrix& operator=(const SparseDOKMatrix &other)
	{
		SparseCSRMatrix tmp(other);
		assign(*this, tmp);
		return *this;
	}

	SparseCSRMatrix& operator=(const SparseIJVMatrix &other)
	{
		SparseCSRMatrix tmp(other);
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

	const MKL_INT* rowPtrs() const
	{
		return &_rowPtrs[0];
	}

	const MKL_INT* columnIndices() const
	{
		return &_columnIndices[0];
	}

	double* values()
	{
		return &_values[0];
	}

	MKL_INT* rowPtrs()
	{
		return &_rowPtrs[0];
	}

	MKL_INT* columnIndices()
	{
		return &_columnIndices[0];
	}

	void dump()
	{
		for(size_t i = 0; i < _values.size(); i++) {
			std::cout << _values[i] << " ";
		}
		std::cout << std::endl;
		for(size_t i = 0; i < _columnIndices.size(); i++) {
			std::cout << _columnIndices[i] << " ";
		}
		std::cout << std::endl;
		for(int i = 0; i <= rows(); i++) {
			std::cout << _rowPtrs[i] << " ";
		}
		std::cout << std::endl;
	}

protected:

	void makeTransposition();

private:

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
