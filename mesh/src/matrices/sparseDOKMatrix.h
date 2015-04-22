#ifndef SPARSEDOKMATRIX_H_
#define SPARSEDOKMATRIX_H_

#include <map>

#include "matrix.h"
#include "denseMatrix.h"
#include "sparseCSRMatrix.h"
#include "sparseIJVMatrix.h"

typedef std::map<MKL_INT, std::map<MKL_INT, double> > MatrixMap;
typedef std::map<MKL_INT, double> ColumnMap;

class DenseMatrix;
class SparseCSRMatrix;
class SparseIJVMatrix;

class SparseDOKMatrix: public EditableMatrix
{

public:

	SparseDOKMatrix(MatrixType type, MKL_INT rowsAndCols): EditableMatrix(type, rowsAndCols, rowsAndCols) {};
	SparseDOKMatrix(MatrixType type, MKL_INT rows, MKL_INT cols): EditableMatrix(type, rows, cols) {};
	SparseDOKMatrix(MKL_INT rows, MKL_INT cols): EditableMatrix(Matrix::GENERAL, rows, cols) {};

	SparseDOKMatrix(const DenseMatrix &other);
	SparseDOKMatrix(const SparseCSRMatrix &other);
	SparseDOKMatrix(const SparseIJVMatrix &other);

	SparseDOKMatrix& operator=(const DenseMatrix &other)
	{
		SparseDOKMatrix tmp(other);
		assign(*this, tmp);
		return *this;
	}

	SparseDOKMatrix& operator=(const SparseCSRMatrix &other)
	{
		SparseDOKMatrix tmp(other);
		assign(*this, tmp);
		return *this;
	}

	SparseDOKMatrix& operator=(const SparseIJVMatrix &other)
	{
		SparseDOKMatrix tmp(other);
		assign(*this, tmp);
		return *this;
	}

	MKL_INT nonZeroValues() const;

	double operator()(MKL_INT row, MKL_INT column) const
	{
		arrange(row, column);

		MatrixMap::const_iterator row_it = _values.find(row);
		if (row_it == _values.end()) {
			return 0;
		}
		ColumnMap::const_iterator column_it = row_it->second.find(column);
		if (column_it == row_it->second.end()) {
			return 0;
		}
		return column_it->second;
	};

	double& operator () (MKL_INT row, MKL_INT column)
	{
		arrange(row, column);
		return _values[row][column];
	}

	MKL_INT removeZeros();

	const std::map<MKL_INT, std::map<MKL_INT, double> >& values() const
	{
		return _values;
	}

	std::map<MKL_INT, std::map<MKL_INT, double> >& values()
	{
		return _values;
	}

	void setRows (MKL_INT rows_new) {
		_rows = rows_new;
	}

	void setCols (MKL_INT cols_new) {
		_cols = cols_new;
	}

protected:

	void makeTransposition();

private:

	static void assign(SparseDOKMatrix &m1, SparseDOKMatrix &m2)
	{
		EditableMatrix::assign(m1, m2);
		m1._values.swap(m2._values);
	}

	MatrixMap _values;
};

#endif /* SPARSEDOKVMATRIX_H_ */
