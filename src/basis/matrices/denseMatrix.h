#ifndef DENSEMATRIX_H_
#define DENSEMATRIX_H_

#include <cstdio>
#include <algorithm>

#include "matrix.h"
template<typename Tindices> class SparseDOKMatrix;
template<typename Tindices> class SparseIJVMatrix;
template<typename Tindices> class SparseCSRMatrix;
template<typename Tindices> class SparseVVPMatrix;


namespace espreso {

#define DenseMatrixIndexing Matrix::ZeroBased

class DenseMatrix: public Matrix
{

public:

	DenseMatrix(): Matrix(DenseMatrixIndexing) {};
	DenseMatrix(size_t rows, size_t columns)
		: Matrix(rows, columns, DenseMatrixIndexing), _values(rows * columns, 0) {};

	template<typename Tindices> DenseMatrix(const SparseDOKMatrix<Tindices> &other);
	template<typename Tindices> DenseMatrix(const SparseCSRMatrix<Tindices> &other);
	template<typename Tindices> DenseMatrix(const SparseIJVMatrix<Tindices> &other);

	template<typename Tindices> DenseMatrix& operator=(const SparseDOKMatrix<Tindices> &other);
	template<typename Tindices> DenseMatrix& operator=(const SparseCSRMatrix<Tindices> &other);
	template<typename Tindices> DenseMatrix& operator=(const SparseIJVMatrix<Tindices> &other);

	DenseMatrix& operator=(double value);
	DenseMatrix operator*(DenseMatrix &M);

	void multiply(
			const DenseMatrix &A, const DenseMatrix &B,
			double alfa = 1, double beta = 0,
			bool transposeA = false, bool transposeB = false);

	void resize(size_t rows, size_t columns);
	void transpose();

	size_t nonZeroValues() const
	{
		return std::count_if(_values.begin(), _values.end(), NonZeroValue());
	}

	double operator()(size_t row, size_t column) const
	{
		return _values[row * _columns + column];
	};

	double& operator()(size_t row, size_t column)
	{
		return _values[row * _columns + column];
	}

	double get(size_t row, size_t column) const
	{
		return _values[row * _columns + column];
	}

	void set(size_t row, size_t column, double value)
	{
		_values[row * _columns + column] = value;
	}

	const double* values() const
	{
		return &_values[0];
	}

	double* values()
	{
		return &_values[0];
	}

	double norm() const
	{
		double n = 0;
		for (size_t i = 0; i < rows(); i++) {
			for (size_t j = 0; j < columns(); j++) {
				n += get(i, j) * get(i, j);
			}
		}
		return sqrt(n);
	}

private:

	static void assign(DenseMatrix &m1, DenseMatrix &m2)
	{
		Matrix::assign(m1, m2);
		m1._values.swap(m2._values);
	}

	std::vector<double> _values;
};

}

#include "denseMatrix.hpp"

#endif /* DENSEMATRIX_H_ */
