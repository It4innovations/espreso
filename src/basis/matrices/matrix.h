#ifndef MATRIX_H_
#define MATRIX_H_

#include <iostream>
#include <cstdio>
#include <vector>
#include <cmath>

#include "mkl_spblas.h"
#include "mkl_blas.h"
#include "mkl_cblas.h"
#include "mkl_lapacke.h"
#include "mkl_trans.h"

#include "../logging/logging.h"

namespace espreso {

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
	enum Indexing {
		ZeroBased,
		OneBased
	};

	size_t rows() const
	{
		return _rows;
	}

	size_t columns() const
	{
		return _columns;
	}

	virtual void resize(size_t rows, size_t columns)
	{
		_rows = rows;
		_columns = columns;
	}

	Indexing indexing() const
	{
		return _indexing;
	}

	virtual double operator()(size_t row, size_t column) const = 0;
	virtual double& operator()(size_t row, size_t column) = 0;

	virtual double get(size_t row, size_t column) const = 0;
	virtual void set(size_t row, size_t column, double value) = 0;

	virtual double norm() const = 0;

	virtual void transpose() = 0;
	virtual size_t nonZeroValues() const = 0;

	virtual ~Matrix() { };

	friend std::ostream& operator<<(std::ostream& os, const Matrix &m);

protected:

	Matrix(Indexing indexing): _rows(0), _columns(0), _indexing(indexing) { };
	Matrix(size_t rows, size_t columns, Indexing indexing)
		: _rows(rows), _columns(columns), _indexing(indexing) { };


	static void assign(Matrix &m1, Matrix &m2)
	{
		m1._rows = m2._rows;
		m1._columns = m2._columns;
	}

	size_t _rows;
	size_t _columns;

	Indexing _indexing;

	static NonZeroValue nonZero;
};

}

#endif /* MATRIX_H_ */
