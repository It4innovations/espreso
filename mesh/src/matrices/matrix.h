#ifndef MATRIX_H_
#define MATRIX_H_

#include <iostream>
#include <cstdio>
#include <vector>
#include "mkl.h"

class Matrix
{

public:

	enum MatrixType {
		GENERAL,
		SYMETRIC
	};

	friend std::ostream& operator<<(std::ostream& os, const Matrix &m);

	virtual double operator()(MKL_INT row, MKL_INT column) const = 0;

	virtual MKL_INT nonZeroValues() const = 0;

	void transpose();

	MKL_INT rows() const
	{
		return _rows;
	}
	MKL_INT columns() const
	{
		return _cols;
	}

	MatrixType type() const
	{
		return _type;
	}

	virtual ~Matrix() { };

protected:

	Matrix(MatrixType type): _type(type), _rows(0), _cols(0) {};
	Matrix(MatrixType type, MKL_INT rows, MKL_INT cols);

	virtual void makeTransposition() = 0;

	void arrange(MKL_INT &row, MKL_INT &column) const
	{
		if (_type == Matrix::SYMETRIC) {
			if (row > column) {
				MKL_INT tmp;
				tmp = row;
				row = column;
				column = tmp;
			}
		}
	}

	static void assign(Matrix &m1, Matrix &m2)
	{
		m1._type = m2._type;
		m1._rows = m2._rows;
		m1._cols = m2._cols;
	}

	MatrixType _type;

	// Variables
	MKL_INT _rows;		// number of rows
	MKL_INT _cols;		// number of columns
};


class EditableMatrix: public Matrix
{

public:

	virtual double operator()(MKL_INT row, MKL_INT column) const = 0;
	virtual double& operator()(MKL_INT row, MKL_INT column) = 0;

	virtual MKL_INT nonZeroValues() const = 0;

	virtual ~EditableMatrix() { };

protected:

	EditableMatrix(MatrixType type, MKL_INT rows, MKL_INT cols) : Matrix(type, rows, cols) { };

	virtual void makeTransposition() = 0;

	static void assign(EditableMatrix &m1, EditableMatrix &m2)
	{
		Matrix::assign(m1, m2);
	}

private:


};

#endif /* MATRIX_H_ */
