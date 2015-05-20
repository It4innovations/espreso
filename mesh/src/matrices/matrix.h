#ifndef MATRIX_H_
#define MATRIX_H_

#include <iostream>
#include <cstdio>
#include <vector>
#include <cmath>
#include "mkl.h"

struct NonZeroValue
{
	bool operator()(double value)
	{
		return fabs(value) > 0;
	}
};

class Matrix
{

public:
	size_t rows() const
	{
		return _rows;
	}

	size_t columns() const
	{
		return _columns;
	}

	void resize(size_t rows, size_t columns)
	{
		_rows = rows;
		_columns = columns;
	}

	virtual double operator()(size_t row, size_t column) const = 0;
	virtual double& operator()(size_t row, size_t column) = 0;

	virtual double get(size_t row, size_t column) const = 0;
	virtual void set(size_t row, size_t column, double value) = 0;

	virtual void transpose() = 0;
	virtual size_t nonZeroValues() const = 0;

	virtual ~Matrix() { };

	friend std::ostream& operator<<(std::ostream& os, const Matrix &m);

protected:

	Matrix(): _rows(0), _columns(0) { };
	Matrix(size_t rows, size_t columns): _rows(rows), _columns(columns) { };



	static void assign(Matrix &m1, Matrix &m2)
	{
		m1._rows = m2._rows;
		m1._columns = m2._columns;
	}

	size_t _rows;
	size_t _columns;

	static NonZeroValue nonZero;
};

#endif /* MATRIX_H_ */
