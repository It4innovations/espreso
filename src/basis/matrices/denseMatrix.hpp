#include "denseMatrix.h"

namespace espreso {

template<typename Tindices>
DenseMatrix::DenseMatrix(const SparseDOKMatrix<Tindices> &other): Matrix(other.rows(), other.columns(), DenseMatrixIndexing)
{
	std::vector<double>(rows() * columns(), 0).swap(_values);

	typename std::map<Tindices, std::map<Tindices, double> >::const_iterator row;
	for(row = other.values().begin(); row != other.values().end(); ++row) {
		typename std::map<Tindices, double>::const_iterator column;
		for(column = row->second.begin(); column != row->second.end(); ++column) {
			set(row->first - other.indexing(), column->first - other.indexing(), column->second);
		}
	}
}

template<typename Tindices>
DenseMatrix::DenseMatrix(const SparseCSRMatrix<Tindices> &other): Matrix(other.rows(), other.columns(), DenseMatrixIndexing)
{
	std::vector<double>(rows() * columns(), 0).swap(_values);

	const Tindices *rowPtrs = other.rowPtrs();
	const Tindices *columnIndices = other.columnIndices();
	const double *values = other.values();

	for (Tindices r = 0; r < _rows; r++) {
		for (Tindices i = rowPtrs[r]; i < rowPtrs[r + 1]; i++) {
			set(r, columnIndices[i - other.indexing()] - other.indexing(), values[i - other.indexing()]);
		}
	}
}

template<typename Tindices>
DenseMatrix::DenseMatrix(const SparseIJVMatrix<Tindices> &other): Matrix(other.rows(), other.columns(), DenseMatrixIndexing)
{
	std::vector<double>(rows() * columns(), 0).swap(_values);

	const Tindices *rowIndices = other.rowIndices();
	const Tindices *columnIndices = other.columnIndices();
	const double *values = other.values();

	for (size_t i = 0; i < other.nonZeroValues(); i++) {
		set(rowIndices[i] - other.indexing(), columnIndices[i] - other.indexing(), values[i]);
	}
}

template<typename Tindices>
DenseMatrix& DenseMatrix::operator=(const SparseDOKMatrix<Tindices> &other)
{
	DenseMatrix tmp(other);
	assign(*this, tmp);
	return *this;
}

template<typename Tindices>
DenseMatrix& DenseMatrix::operator=(const SparseCSRMatrix<Tindices> &other)
{
	DenseMatrix tmp(other);
	assign(*this, tmp);
	return *this;
}

template<typename Tindices>
DenseMatrix& DenseMatrix::operator=(const SparseIJVMatrix<Tindices> &other)
{
	DenseMatrix tmp(other);
	assign(*this, tmp);
	return *this;
}

}
