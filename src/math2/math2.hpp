
#ifndef SRC_MATH2_MATH2_HPP_
#define SRC_MATH2_MATH2_HPP_

#include "math2.h"

namespace espreso {
namespace math {

template <typename T> void copy(Matrix_Dense<T>  &x, const Matrix_CSR<T>   &y)
{
	eslog::error("call empty function copy\n");
}

template <typename T> void copy(Matrix_Dense<T>  &x, const Matrix_IJV<T>   &y)
{
	eslog::error("call empty function copy\n");
}

template <typename T> void copy(Matrix_CSR<T>    &x, const Matrix_Dense<T> &y)
{
	eslog::error("call empty function copy\n");
}

template <typename T> void copy(Matrix_CSR<T>    &x, const Matrix_IJV<T>   &y)
{
	eslog::error("call empty function copy\n");
}

template <typename T> void copy(Matrix_IJV<T>    &x, const Matrix_Dense<T> &y)
{
	eslog::error("call empty function copy\n");
}

template <typename T> void copy(Matrix_IJV<T>    &x, const Matrix_CSR<T>   &y)
{
	eslog::error("call empty function copy\n");
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
	if (x.nnz < y.nnz) {
		for (esint i = 0; i < x.nnz / size; ++i) {
			for (esint j = 0; j < size; ++j) {
				x.vals[i * size + j] = y.vals[i * step + offset + j];
			}
		}
	} else {
		for (esint i = 0; i < y.nnz / size; ++i) {
			for (esint j = 0; j < size; ++j) {
				x.vals[i * step + offset + j] = y.vals[i * size + j];
			}
		}
	}
}

template <typename T> void copy(Vector_Dense<T>  &x, const Vector_Sparse<T> &y, int offset, int size, int step)
{
	eslog::error("call empty function copy\n");
}

template <typename T> void copy(Vector_Sparse<T> &x, const Vector_Dense<T>  &y, int offset, int size, int step)
{
	eslog::error("call empty function copy\n");
}

template <typename T> void copy(Matrix_Dense<T> &x, const Matrix_Dense<T> &y, int rowOffset, int colOffset, int size, int step)
{
	eslog::error("call empty function copy\n");
}

template <typename T> void copy(Matrix_Dense<T> &x, const Matrix_CSR<T> &y, int rowOffset, int colOffset, int size, int step)
{
	eslog::error("call empty function copy\n");
}

template <typename T> void copy(Matrix_Dense<T> &x, const Matrix_IJV<T> &y, int rowOffset, int colOffset, int size, int step)
{
	eslog::error("call empty function copy\n");
}

template <typename T> void copy(Matrix_CSR<T> &x, const Matrix_Dense<T> &y, int rowOffset, int colOffset, int size, int step)
{
	eslog::error("call empty function copy\n");
}

template <typename T> void copy(Matrix_CSR<T> &x, const Matrix_CSR<T> &y, int rowOffset, int colOffset, int size, int step)
{
	if (x.nnz < y.nnz) {
		eslog::error("call empty function copy\n");
//		for (esint r = 0; r < y.nrows; r++) {
//			for (esint c = 0; c < y.rows[r + 1] - y.rows[r]; c++) {
//				esint tr = step * (r / size) + rowOffset + r % size;
//				esint tc = step * (c / size) + colOffset + c % size;
//				x.vals[x.rows[tr] + tc - Indexing::CSR] = y.vals[y.rows[r] + c - Indexing::CSR];
//			}
//		}
	} else {
		for (esint r = 0; r < y.nrows; r++) {
			for (esint c = 0; c < y.rows[r + 1] - y.rows[r]; c++) {
				esint tr = step * (r / size) + rowOffset + r % size;
				esint tc = step * (c / size) + colOffset + c % size;
				x.vals[x.rows[tr] + tc - Indexing::CSR] = y.vals[y.rows[r] + c - Indexing::CSR];
			}
		}
	}
}

template <typename T> void copy(Matrix_CSR<T> &x, const Matrix_IJV<T> &y, int rowOffset, int colOffset, int size, int step)
{
	eslog::error("call empty function copy\n");
}

template <typename T> void copy(Matrix_IJV<T> &x, const Matrix_Dense<T> &y, int rowOffset, int colOffset, int size, int step)
{
	eslog::error("call empty function copy\n");
}

template <typename T> void copy(Matrix_IJV<T> &x, const Matrix_CSR<T> &y, int rowOffset, int colOffset, int size, int step)
{
	eslog::error("call empty function copy\n");
}

template <typename T> void copy(Matrix_IJV<T> &x, const Matrix_IJV<T> &y, int rowOffset, int colOffset, int size, int step)
{
	eslog::error("call empty function copy\n");
}

template <typename T> void copy(Matrix_Dense<std::complex<T> > &x, const int offsetX, const Matrix_CSR<T>   &y)
{
	eslog::error("call empty function copy\n");
}

template <typename T> void copy(Matrix_Dense<std::complex<T> > &x, const int offsetX, const Matrix_IJV<T>   &y)
{
	eslog::error("call empty function copy\n");
}

template <typename T> void copy(Matrix_CSR<std::complex<T> >   &x, const int offsetX, const Matrix_Dense<T> &y)
{
	eslog::error("call empty function copy\n");
}

template <typename T> void copy(Matrix_CSR<std::complex<T> >   &x, const int offsetX, const Matrix_IJV<T>   &y)
{
	eslog::error("call empty function copy\n");
}

template <typename T> void copy(Matrix_IJV<std::complex<T> >   &x, const int offsetX, const Matrix_Dense<T> &y)
{
	eslog::error("call empty function copy\n");
}

template <typename T> void copy(Matrix_IJV<std::complex<T> >   &x, const int offsetX, const Matrix_CSR<T>   &y)
{
	eslog::error("call empty function copy\n");
}

template <typename T> void add(Matrix_Dense<T>  &x, const T &alpha, const Matrix_CSR<T>   &y)
{
	eslog::error("call empty function add\n");
}

template <typename T> void add(Matrix_Dense<T>  &x, const T &alpha, const Matrix_IJV<T>   &y)
{
	eslog::error("call empty function add\n");
}

template <typename T> void add(Matrix_CSR<T>    &x, const T &alpha, const Matrix_Dense<T> &y)
{
	eslog::error("call empty function add\n");
}

template <typename T> void add(Matrix_CSR<T>    &x, const T &alpha, const Matrix_IJV<T>   &y)
{
	eslog::error("call empty function add\n");
}

template <typename T> void add(Matrix_IJV<T>    &x, const T &alpha, const Matrix_Dense<T> &y)
{
	eslog::error("call empty function add\n");
}

template <typename T> void add(Matrix_IJV<T>    &x, const T &alpha, const Matrix_CSR<T>   &y)
{
	eslog::error("call empty function add\n");
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

template <typename T> void add(Matrix_Dense<std::complex<T>> &x, const int offsetX, const T &alpha, const Matrix_CSR<T>   &y)
{
	eslog::error("call empty function add\n");
}

template <typename T> void add(Matrix_Dense<std::complex<T>> &x, const int offsetX, const T &alpha, const Matrix_IJV<T>   &y)
{
	eslog::error("call empty function add\n");
}

template <typename T> void add(Matrix_CSR<std::complex<T>>   &x, const int offsetX, const T &alpha, const Matrix_Dense<T> &y)
{
	eslog::error("call empty function add\n");
}

template <typename T> void add(Matrix_CSR<std::complex<T>>   &x, const int offsetX, const T &alpha, const Matrix_IJV<T>   &y)
{
	eslog::error("call empty function add\n");
}

template <typename T> void add(Matrix_IJV<std::complex<T>>   &x, const int offsetX, const T &alpha, const Matrix_Dense<T> &y)
{
	eslog::error("call empty function add\n");
}

template <typename T> void add(Matrix_IJV<std::complex<T>>   &x, const int offsetX, const T &alpha, const Matrix_CSR<T>   &y)
{
	eslog::error("call empty function add\n");
}

template <typename T> void add(Matrix_Dense<T> &x, const T &beta, const Matrix_Dense<T> &y, int rowOffset, int colOffset, int size, int step)
{
	eslog::error("call empty function add\n");
}

template <typename T> void add(Matrix_Dense<T> &x, const T &beta, const Matrix_CSR<T> &y, int rowOffset, int colOffset, int size, int step)
{
	eslog::error("call empty function add\n");
}

template <typename T> void add(Matrix_Dense<T> &x, const T &beta, const Matrix_IJV<T> &y, int rowOffset, int colOffset, int size, int step)
{
	eslog::error("call empty function add\n");
}

template <typename T> void add(Matrix_CSR<T> &x, const T &beta, const Matrix_Dense<T> &y, int rowOffset, int colOffset, int size, int step)
{
	eslog::error("call empty function add\n");
}

template <typename T> void add(Matrix_CSR<T> &x, const T &beta, const Matrix_CSR<T> &y, int rowOffset, int colOffset, int size, int step)
{
	// assume the same pattern for y and z
	// x has 2 times larger pattern
	for (esint r = 0; r < y.nrows; r++) {
		for (esint c = 0; c < y.rows[r + 1] - y.rows[r]; c++) {
			esint tr = step * (r / size) + rowOffset + r % size;
			esint tc = step * (c / size) + colOffset + c % size;
			x.vals[x.rows[tr] + tc - Indexing::CSR] += beta * y.vals[y.rows[r] + c - Indexing::CSR];
		}
	}
}

template <typename T> void add(Matrix_CSR<T> &x, const T &beta, const Matrix_IJV<T> &y, int rowOffset, int colOffset, int size, int step)
{
	eslog::error("call empty function add\n");
}

template <typename T> void add(Matrix_IJV<T> &x, const T &beta, const Matrix_Dense<T> &y, int rowOffset, int colOffset, int size, int step)
{
	eslog::error("call empty function add\n");
}

template <typename T> void add(Matrix_IJV<T> &x, const T &beta, const Matrix_CSR<T> &y, int rowOffset, int colOffset, int size, int step)
{
	eslog::error("call empty function add\n");
}

template <typename T> void add(Matrix_IJV<T> &x, const T &beta, const Matrix_IJV<T> &y, int rowOffset, int colOffset, int size, int step)
{
	eslog::error("call empty function add\n");
}

template <typename T> void orthonormalize(Matrix_Dense<T> &m)
{
	for (esint r = 0; r < m.nrows; ++r) {
		for (esint rr = 0; rr < r; ++rr) {
			double scale = math::dot(m.ncols, m.vals + rr * m.ncols, 1, m.vals + r * m.ncols, 1) / math::dot(m.ncols, m.vals + rr * m.ncols, 1, m.vals + rr * m.ncols, 1);
			math::add(m.ncols, m.vals + r * m.ncols, 1, -scale, m.vals + rr * m.ncols, 1);
		}
		math::scale(m.ncols, 1. / math::norm(m.ncols, m.vals + r * m.ncols, 1), m.vals + r * m.ncols, 1);
	}
}

}
}

#endif /* SRC_MATH2_MATH2_HPP_ */
