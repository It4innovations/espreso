
#ifndef SRC_ANALYSIS_MATH_MATH_PHYSICS_COPY_H_
#define SRC_ANALYSIS_MATH_MATH_PHYSICS_COPY_H_

#include "analysis/math/math.physics.h"

namespace espreso {
namespace math {


template <typename T> void copy(Matrix_Dense<T>  &x, const Matrix_CSR<T>   &y)
{
	if (x.nrows != y.nrows || x.ncols != y.ncols) {
		eslog::error("matrices have incorrect dimensions.\n");
	}
	math::set(x, T{0});
	for (esint r = 0; r < y.nrows; ++r) {
		for (esint c = y.rows[r]; c < y.rows[r + 1]; ++c) {
			x.vals[r * x.ncols + y.cols[c]] = y.vals[c];
		}
	}
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

}
}

#endif /* SRC_ANALYSIS_MATH_MATH_PHYSICS_COPY_H_ */
