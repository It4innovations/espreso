#ifndef SPARSEIJVMATRIX_H_
#define SPARSEIJVMATRIX_H_

#include "matrix.h"
#include "denseMatrix.h"
#include "sparseDOKMatrix.h"
#include "sparseCSRMatrix.h"
#include "sparseVVPMatrix.h"

namespace espreso {

class DenseMatrix;
template<typename Tindices> class SparseDOKMatrix;
template<typename Tindices> class SparseCSRMatrix;
template<typename Tindices> class SparseVVPMatrix;

#define IJVMatrixIndexing Matrix::OneBased

template<typename Tindices>
class SparseIJVMatrix: public Matrix
{

public:

	template<typename TIndices>
	friend std::ostream& operator<<(std::ostream& os, const SparseIJVMatrix<TIndices> &m);

	SparseIJVMatrix(): Matrix(IJVMatrixIndexing) { };
	SparseIJVMatrix(size_t rows, size_t columns): Matrix(rows, columns, IJVMatrixIndexing) { };

	SparseIJVMatrix(const DenseMatrix &other);
	SparseIJVMatrix(const SparseDOKMatrix<Tindices> &other);
	SparseIJVMatrix(const SparseCSRMatrix<Tindices> &other);
	SparseIJVMatrix(SparseVVPMatrix<Tindices> &other);

	SparseIJVMatrix<Tindices>& operator=(const DenseMatrix &other);
	SparseIJVMatrix<Tindices>& operator=(const SparseDOKMatrix<Tindices> &other);
	SparseIJVMatrix<Tindices>& operator=(const SparseCSRMatrix<Tindices> &other);
	SparseIJVMatrix<Tindices>& operator=(SparseVVPMatrix<Tindices> &other);

	void reserve(size_t size);
	void transpose();

	void sort()
	{
		// last element is part of the array
		sort(0, _values.size() - 1);
	}

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
		for(eslocal i = 0; i < _rowIndices.size(); i++)
		{
			if (_rowIndices[i] == row + _indexing && _columnIndices[i] == column + _indexing) {
				return _values[i];
			}
		}
		return 0;
	}

	std::vector<Tindices>& rowIndices()
	{
		return _rowIndices;
	}

	const std::vector<Tindices>& rowIndices() const
	{
		return _rowIndices;
	}

	std::vector<Tindices>& columnIndices()
	{
		return _columnIndices;
	}

	const std::vector<Tindices>& columnIndices() const
	{
		return _columnIndices;
	}

	std::vector<double>& values()
	{
		return _values;
	}

	const std::vector<double>& values() const
	{
		return _values;
	}

	void ShiftRowIndex(eslocal offset) {
		for (int i = 0; i < _rowIndices.size(); i++) {
			_rowIndices[i] += offset;
		}

		if (_rowIndices.size() && _rowIndices.back() > _rows) {
			_rows = _rowIndices.back();
		}
	}

	void AppendMatrix(SparseIJVMatrix &inputMatrix) {
		_rowIndices.insert(_rowIndices.end(), inputMatrix._rowIndices.begin(), inputMatrix._rowIndices.end());
		_columnIndices.insert(_columnIndices.end(), inputMatrix._columnIndices.begin(), inputMatrix._columnIndices.end());
		_values.insert(_values.end(), inputMatrix._values.begin(), inputMatrix._values.end());

		_rows    = _rows > inputMatrix.rows() ? _rows : inputMatrix.rows();
		_columns = _columns > inputMatrix.columns() ? _columns : inputMatrix.columns();

		//sort();
	}

	void push(size_t row, size_t column, double value)
	{
		_rowIndices.push_back(row + _indexing);
		_columnIndices.push_back(column + _indexing);
		_values.push_back(value);
	}

	double& operator()(size_t row, size_t column)
	{
		_rowIndices.push_back(row + _indexing);
		_columnIndices.push_back(column + _indexing);
		_values.push_back(0);
		return _values.back();
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

	void sort(size_t begin, size_t end);

	static void assign(SparseIJVMatrix<Tindices> &m1, SparseIJVMatrix<Tindices> &m2)
	{
		Matrix::assign(m1, m2);
		m1._rowIndices.swap(m2._rowIndices);
		m1._columnIndices.swap(m2._columnIndices);
		m1._values.swap(m2._values);
	}

	void set(size_t row, size_t column, double value)
	{
		if (Matrix::nonZero(value)) {
			for(size_t i = 0; i < _rowIndices.size(); i++)
			{
				if (_rowIndices[i] == row + _indexing && _columnIndices[i] == column + _indexing) {
					_values[i] = value;
				}
			}
			_rowIndices.push_back(row + _indexing);
			_columnIndices.push_back(column + _indexing);
			_values.push_back(value);
		}
	}

	inline Tindices compare(size_t i, size_t j)
	{
		if (_rowIndices[i] == _rowIndices[j]) {
			return _columnIndices[i] - _columnIndices[j];
		} else {
			return _rowIndices[i] - _rowIndices[j];
		}
	}

	inline void swap(size_t i, size_t j)
	{
		Tindices t = _columnIndices[i];
		_columnIndices[i] = _columnIndices[j];
		_columnIndices[j] = t;

		t = _rowIndices[i];
		_rowIndices[i] = _rowIndices[j];
		_rowIndices[j] = t;

		double td = _values[i];
		_values[i] = _values[j];
		_values[j] = td;
	}

	// Sparse COO data
	std::vector<Tindices> _rowIndices;
	std::vector<Tindices> _columnIndices;
	std::vector<double> _values;
};

}

#include "sparseIJVMatrix.hpp"

#endif /* SPARSEIJVMATRIX_H_ */
