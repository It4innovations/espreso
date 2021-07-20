
#ifndef SRC_MATH2_MATH2_HPP_
#define SRC_MATH2_MATH2_HPP_

#include "math2.h"

namespace espreso {
namespace math {

template <typename T> void multiplyPattern(Vector_Dense<T> &x, const Vector_Dense<T> &y, int size, int multipicity)
{

}

template <typename T> void multiplyPattern(Vector_Sparse<T> &x, const Vector_Sparse<T> &y, int size, int multipicity)
{
	x.resize(multipicity * y.size, multipicity * y.nnz);
	for (esint i = 0; i < y.nnz / size; ++i) {
		for (int j = 0; j < size; ++j) {
			for (int m = 0; m < multipicity; ++m) {
				x.indices[i * size * multipicity + j * multipicity + m] = multipicity * y.indices[i * size + j] + m;
			}
		}
	}
}

template <typename T> void multiplyPattern(Matrix_Dense<T> &x, const Matrix_Dense<T> &y, int size, int multipicity)
{

}

template <typename T> void multiplyPattern(Matrix_CSR<T> &x, const Matrix_CSR<T> &y, int size, int multipicity)
{

}

template <typename T> void multiplyPattern(Matrix_IJV<T> &x, const Matrix_IJV<T> &y, int size, int multipicity)
{

}

template <typename T> void copy(Vector_Dense<T> &x, const Vector_Dense<T> &y, int offset, int size, int step)
{
	if (x.size < y.size) {
		for (esint i = 0; i < x.size / size; ++i) {
			for (esint j = 0; j < size; ++j) {
				x.vals[i * size + j] = y.vals[i * step + offset + j];
			}
		}
	} else {
		for (esint i = 0; i < y.size / size; ++i) {
			for (esint j = 0; j < size; ++j) {
				x.vals[i * step + offset + j] = y.vals[i * size + j];
			}
		}
	}
}

template <typename T> void copy(Vector_Sparse<T> &x, const Vector_Sparse<T> &y, int offset, int size, int step)
{

}

template <typename T> void copy(Matrix_Dense<T> &x, const Matrix_Dense<T> &y, int rowOffset, int colOffset, int size, int step)
{

}

template <typename T> void copy(Matrix_CSR<T> &x, const Matrix_CSR<T> &y, int rowOffset, int colOffset, int size, int step)
{

}

template <typename T> void copy(Matrix_IJV<T> &x, const Matrix_IJV<T> &y, int rowOffset, int colOffset, int size, int step)
{

}

template <typename T> void add(Vector_Dense<T> &x, const T &beta, const Vector_Dense<T> &y, int offset, int size, int step)
{
	for (esint i = 0; i < y.size / size; ++i) {
		for (esint j = 0; j < size; ++j) {
			x.vals[i * step + offset + j] += beta * y.vals[i * size + j];
		}
	}
}

template <typename T> void add(Vector_Sparse<T> &x, const T &beta, const Vector_Sparse<T> &y, int offset, int size, int step)
{
	for (esint i = 0; i < y.nnz / size; ++i) {
		for (esint j = 0; j < size; ++j) {
			x.vals[i * step + offset + j] += beta * y.vals[i * size + j];
		}
	}
}

template <typename T> void add(Matrix_Dense<T> &x, const T &beta, const Matrix_Dense<T> &y, int rowOffset, int colOffset, int size, int step)
{

}

template <typename T> void add(Matrix_CSR<T> &x, const T &beta, const Matrix_CSR<T> &y, int rowOffset, int colOffset, int size, int step)
{

}

template <typename T> void add(Matrix_IJV<T> &x, const T &beta, const Matrix_IJV<T> &y, int rowOffset, int colOffset, int size, int step)
{

}

template <typename T> void sum(Vector_Dense<T> &x, const T &alpha, const Vector_Dense<T> &y, const T &beta, const Vector_Dense<T> &z, int offset, int size, int step)
{

}

template <typename T> void sum(Vector_Sparse<T> &x, const T &alpha, const Vector_Sparse<T> &y, const T &beta, const Vector_Sparse<T> &z, int offset, int size, int step)
{

}

template <typename T> void sum(Matrix_Dense<T> &x, const T &alpha, const Matrix_Dense<T> &y, const T &beta, const Matrix_Dense<T> &z, int rowOffset, int colOffset, int size, int step)
{

}

template <typename T> void sum(Matrix_CSR<T> &x, const T &alpha, const Matrix_CSR<T> &y, const T &beta, const Matrix_CSR<T> &z, int rowOffset, int colOffset, int size, int step)
{
	// assume the same pattern for y and z
	// x has 2 times larger pattern
	for (esint r = 0; r < y.nrows; r++) {
	for (esint c = 0; c < y.rows[r + 1] - y.rows[r]; c++) {
		esint tr = step * (r / size) + rowOffset + r % size;
		esint tc = step * (c / size) + colOffset + c % size;
		x.vals[x.rows[tr] + tc - _Matrix_CSR_Pattern::Indexing] = alpha * y.vals[y.rows[r] + c - _Matrix_CSR_Pattern::Indexing] + beta * z.vals[y.rows[r] + c - _Matrix_CSR_Pattern::Indexing];
	}
}
}

template <typename T> void sum(Matrix_IJV<T> &x   , const T &alpha, const Matrix_IJV<T> &y   , const T &beta, const Matrix_IJV<T> &z   , int rowOffset, int colOffset, int size, int step)
{

}

}
}

#endif /* SRC_MATH2_MATH2_HPP_ */
