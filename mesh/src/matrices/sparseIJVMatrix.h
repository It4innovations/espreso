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

#define IJVMatrixIndexing Matrix::OneBased

class SparseIJVMatrix: public Matrix
{

public:

	friend std::ostream& operator<<(std::ostream& os, const SparseIJVMatrix &m);

	SparseIJVMatrix(): Matrix(IJVMatrixIndexing) { };
	SparseIJVMatrix(eslocal rows, eslocal columns): Matrix(rows, columns, IJVMatrixIndexing) { };

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
		for(eslocal i = 0; i < _rowIndices.size(); i++)
		{
			if (_rowIndices[i] == row + _indexing && _columnIndices[i] == column + _indexing) {
				return _values[i];
			}
		}
		return 0;
	}

	const eslocal* rowIndices() const
	{
		return &_rowIndices[0];
	}

	eslocal* rowIndices()
	{
		return &_rowIndices[0];
	}

	std::vector < eslocal > & rowIndices_vec () {
		return _rowIndices;
	}

	const eslocal* columnIndices() const
	{
		return &_columnIndices[0];
	}

	eslocal* columnIndices()
	{
		return &_columnIndices[0];
	}

	std::vector < eslocal > & columnIndices_vec () {
		return _columnIndices;
	}

	const double* values() const
	{
		return &_values[0];
	}

	double* values()
	{
		return &_values[0];
	}

	std::vector < double > & values_vec () {
		return _values;
	}

	void ShiftRowIndex ( eslocal offset ) {
		for (int i = 0; i < _rowIndices.size(); i++)
			_rowIndices[i]+=offset;

		if ( _rowIndices.size() > 0 && _rowIndices[_rowIndices.size()-1] > _rows )
				_rows = _rowIndices[_rowIndices.size()-1];

	}

	void AppendMatrix ( SparseDOKMatrix &inputMatrix) {
		SparseIJVMatrix inpMatIJV;
		inpMatIJV = inputMatrix;

		_rowIndices.insert   (_rowIndices.end(),    inpMatIJV.rowIndices_vec().begin(),    inpMatIJV.rowIndices_vec().end()    );
		_columnIndices.insert(_columnIndices.end(), inpMatIJV.columnIndices_vec().begin(), inpMatIJV.columnIndices_vec().end() );
		_values.insert(       _values.end(),        inpMatIJV.values_vec().begin(),        inpMatIJV.values_vec().end()        );

		_rows    = ( (_rows    > inpMatIJV.rows())    ? _rows    : inpMatIJV.rows()    );
		_columns = ( (_columns > inpMatIJV.columns()) ? _columns : inpMatIJV.columns() );


//		if ( (_rowIndices[_rowIndices.size()-1] + _indexing ) >= _rows ) {
//			_rows = _rowIndices[_rowIndices.size()-1] + _indexing;
//		}
//
//
//		if ( (_columnIndices[_columnIndices.size()-1] +_indexing ) >= _columns ) {
//			_columns = _columnIndices[_columnIndices.size()-1] + _indexing;
//		}

		//TODO:Call IJVSort()

	}



private:

	double& operator()(size_t row, size_t column)
	{
		for(size_t i = 0; i < _rowIndices.size(); i++)
		{
			if (_rowIndices[i] == row + _indexing && _columnIndices[i] == column + _indexing) {
				return _values[i];
			}
		}
		_rowIndices.push_back(row + _indexing);
		_columnIndices.push_back(column + _indexing);
		_values.push_back(0);
		return _values.back();
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

	static void assign(SparseIJVMatrix &m1, SparseIJVMatrix &m2)
	{
		Matrix::assign(m1, m2);
		m1._rowIndices.swap(m2._rowIndices);
		m1._columnIndices.swap(m2._columnIndices);
		m1._values.swap(m2._values);

	}

	// Sparse COO data
	std::vector<eslocal> _rowIndices;
	std::vector<eslocal> _columnIndices;
	std::vector<double> _values;
};

#endif /* SPARSEIJVMATRIX_H_ */
