
#ifndef SRC_MATH2_MATH2_H_
#define SRC_MATH2_MATH2_H_

#include "esinfo/eslog.h"
#include "primitives/vector_dense.h"
#include "primitives/vector_sparse.h"
#include "primitives/matrix_dense.h"
#include "primitives/matrix_csr.h"
#include "primitives/matrix_ijv.h"

#include <complex>

namespace espreso {
namespace math { // interface to wrappers

	// utility functions allowing the Intel inspector-executor model
	template <typename T> void commit(Matrix_Dense<T> &x);
	template <typename T> void commit(Matrix_CSR<T> &x);
	template <typename T> void commit(Matrix_IJV<T> &x);

	template <typename T> void free(Matrix_Dense<T> &x);
	template <typename T> void free(Matrix_CSR<T> &x);
	template <typename T> void free(Matrix_IJV<T> &x);

	template <typename T> void symbolicFactorization(const Matrix_CSR<T> &x);
	template <typename T> void numericalFactorization(const Matrix_CSR<T> &x);
	template <typename T> void solve(const Matrix_CSR<T> &x, Vector_Dense<T> &rhs, Vector_Dense<T> &solution);
	template <typename T> void solve(const Matrix_CSR<T> &x, Matrix_Dense<T> &rhs, Matrix_Dense<T> &solution);

	// y = alpha * A * x + beta * y
	template <typename T> void apply(Vector_Dense<T> &y, const T &alpha, const Matrix_Dense<T> &a, const T &beta, const Vector_Dense<T> &x);
	template <typename T> void apply(Vector_Dense<T> &y, const T &alpha, const Matrix_CSR<T>   &a, const T &beta, const Vector_Dense<T> &x);
	template <typename T> void apply(Vector_Dense<T> &y, const T &alpha, const Matrix_IJV<T>   &a, const T &beta, const Vector_Dense<T> &x);

	// x = value
	template <typename T>
	void set(const esint size, T *x, const int incX, const T &value)
	{
		for (esint i = 0; i < size; i += incX) {
			x[i] = value;
		}
	}

	// x = y
	template <typename T>
	void copy(const esint size, T *x, const int incX, const T *y, const int incY);

	// x *= alpha
	template <typename T>
	void scale(const esint size, const T &alpha, T *x, const int incX);

	// x += alpha * y
	template <typename T>
	void add(const esint size, T *x, const int incX, const T &alpha, const T *y, const int incY);

	template <typename T>
	T dot(const esint size, const T *x, const int incX, const T *y, const int incY);

	template <typename T>
	T norm(const esint size, const T *x, const int incX);

	// x = alpha * y + beta * z
//	template <typename T>
//	void sum(const esint size, T *x, const esint incX, const T &alpha, const T *y, const esint incY);

} // math (interface to wrappers)
} // espreso

namespace espreso {
namespace math {

	template <typename T> void set(Vector_Dense<T>  &x, const T &value) { set(x.size           , x.vals, 1, value); }
	template <typename T> void set(Vector_Sparse<T> &x, const T &value) { set(x.nnz            , x.vals, 1, value); }
	template <typename T> void set(Matrix_Dense<T>  &x, const T &value) { set(x.nrows * x.ncols, x.vals, 1, value); }
	template <typename T> void set(Matrix_CSR<T>    &x, const T &value) { set(x.nnz            , x.vals, 1, value); }
	template <typename T> void set(Matrix_IJV<T>    &x, const T &value) { set(x.nnz            , x.vals, 1, value); }

	template <typename T> void copy(Vector_Dense<T>  &x, const Vector_Dense<T>  &y) { copy(x.size           , x.vals, 1, y.vals, 1); }
	template <typename T> void copy(Vector_Sparse<T> &x, const Vector_Sparse<T> &y) { copy(x.nnz            , x.vals, 1, y.vals, 1); }
	template <typename T> void copy(Matrix_Dense<T>  &x, const Matrix_Dense<T>  &y) { copy(x.nrows * x.ncols, x.vals, 1, y.vals, 1); }
	template <typename T> void copy(Matrix_CSR<T>    &x, const Matrix_CSR<T>    &y) { copy(x.nnz            , x.vals, 1, y.vals, 1); }
	template <typename T> void copy(Matrix_IJV<T>    &x, const Matrix_IJV<T>    &y) { copy(x.nnz            , x.vals, 1, y.vals, 1); }

	template <typename T> void scale(const T &alpha, Vector_Dense<T>  &x) { scale(x.size           , alpha, x.vals, 1); }
	template <typename T> void scale(const T &alpha, Vector_Sparse<T> &x) { scale(x.nnz            , alpha, x.vals, 1); }
	template <typename T> void scale(const T &alpha, Matrix_Dense<T>  &x) { scale(x.nrows * x.ncols, alpha, x.vals, 1); }
	template <typename T> void scale(const T &alpha, Matrix_CSR<T>    &x) { scale(x.nnz            , alpha, x.vals, 1); }
	template <typename T> void scale(const T &alpha, Matrix_IJV<T>    &x) { scale(x.nnz            , alpha, x.vals, 1); }

	template <typename T> void add(Vector_Dense<T>  &x, const T &alpha, const Vector_Dense<T>  &y) { add(x.size           , x.vals, 1, alpha, y.vals, 1); }
	template <typename T> void add(Vector_Sparse<T> &x, const T &alpha, const Vector_Sparse<T> &y) { add(x.nnz            , x.vals, 1, alpha, y.vals, 1); }
	template <typename T> void add(Matrix_Dense<T>  &x, const T &alpha, const Matrix_Dense<T>  &y) { add(x.nrows * x.ncols, x.vals, 1, alpha, y.vals, 1); }
	template <typename T> void add(Matrix_CSR<T>    &x, const T &alpha, const Matrix_CSR<T>    &y) { add(x.nnz            , x.vals, 1, alpha, y.vals, 1); }
	template <typename T> void add(Matrix_IJV<T>    &x, const T &alpha, const Matrix_IJV<T>    &y) { add(x.nnz            , x.vals, 1, alpha, y.vals, 1); }

	template <typename T> void combine(Matrix_CSR<T> &C, const Matrix_CSR<T> &A, const Matrix_CSR<T> &B)
	{
		if (A.nrows != B.nrows || A.ncols != B.ncols) {
			eslog::error("invalid matrices sizes.\n");
		}
		esint nnz = 0;
		for (esint r = 0; r < A.nrows; ++r) {
			esint *beginA = A.cols + A.rows[r    ] - Indexing::CSR;
			esint *endA   = A.cols + A.rows[r + 1] - Indexing::CSR;
			esint *beginB = B.cols + B.rows[r    ] - Indexing::CSR;
			esint *endB   = B.cols + B.rows[r + 1] - Indexing::CSR;
			while (true) {
				if (beginA == endA) { nnz += endB - beginB; break; }
				if (beginB == endB) { nnz += endA - beginA; break; }
				if (*beginA == *beginB) {
					++beginA; ++beginB;
				} else {
					if (*beginA < *beginB) {
						++beginA;
					} else {
						++beginB;
					}
				}
				++nnz;
			}
		}
		C.resize(A.nrows, A.ncols, nnz);
		esint *r = C.rows, *c = C.cols;
		for (esint i = 0; i < A.nrows; ++i, ++r) {
			*r = c - C.cols + Indexing::CSR;
			esint *bA = A.cols + A.rows[i    ] - Indexing::CSR;
			esint *eA = A.cols + A.rows[i + 1] - Indexing::CSR;
			esint *bB = B.cols + B.rows[i    ] - Indexing::CSR;
			esint *eB = B.cols + B.rows[i + 1] - Indexing::CSR;
			while (true) {
				if (bA == eA) { while (bB != eB) { *c++ = *bB++;} break; }
				if (bB == eB) { while (bA != eA) { *c++ = *bA++;} break; }
				if (*bA == *bB) {
					*c++ = *bA++; bB++;
				} else {
					if (*bA < *bB) {
						*c++ = *bA++;
					} else {
						*c++ = *bB++;
					}
				}
			}
		}
		*r = c - C.cols + Indexing::CSR;
	}

	template <typename T> void sumCombined(Matrix_CSR<T> &C, const T &alpha, const Matrix_CSR<T> &A, const Matrix_CSR<T> &B)
	{
		if (A.nrows != B.nrows || A.ncols != B.ncols) {
			eslog::error("invalid matrices sizes.\n");
		}
		for (esint r = 0; r < A.nrows; ++r) {
			esint *bA = A.cols + A.rows[r    ] - Indexing::CSR;
			esint *eA = A.cols + A.rows[r + 1] - Indexing::CSR;
			esint *bB = B.cols + B.rows[r    ] - Indexing::CSR;
			esint *eB = B.cols + B.rows[r + 1] - Indexing::CSR;
			esint *bC = C.cols + C.rows[r    ] - Indexing::CSR;
			esint *eC = C.cols + C.rows[r + 1] - Indexing::CSR;
			while (bC != eC) {
				C.vals[bC - C.cols] = 0;
				if (bA != eA && *bC == *bA) { C.vals[bC - C.cols] += A.vals[bA++ - A.cols]; }
				if (bB != eB && *bC == *bB) { C.vals[bC - C.cols] += B.vals[bB++ - B.cols]; }
				++bC;
			}
		}
	}

	template <typename T> T dot(const Vector_Dense<T>  &x, const Vector_Dense<T>  &y) { return dot(x.size, x.vals, 1, y.vals, 1); }
	template <typename T> T dot(const Vector_Sparse<T> &x, const Vector_Sparse<T> &y) { return dot(x.nnz , x.vals, 1, y.vals, 1); }

	template <typename T> T norm(const Vector_Dense<T>  &x) { return norm(x.size, x.vals, 1); }
	template <typename T> T norm(const Vector_Sparse<T> &x) { return norm(x.nnz , x.vals, 1); }

	template <typename T> T max(const Vector_Dense<T> &x)
	{
		T max = x.vals[0];
		for (esint i = 0; i < x.size; ++i) {
			if (max < x.vals[i]) {
				max = x.vals[i];
			}
		}
		return max;
	}

	template <typename T> T getDiagonalMax(const Matrix_CSR<T> &m)
	{
		T max = m.vals[0];
		for (esint r = 0; r < m.nrows; ++r) {
			esint *c = m.cols + m.rows[r] - Indexing::CSR;
			while (*c != r + Indexing::CSR) { ++c; }
			max = std::max(max, m.vals[c - m.cols]);
		}
		return max;
	}

	template <typename T> void orthonormalize(Matrix_Dense<T> &m);

	template <class T> void store(const T &x, const char* file);

} // math
} // espreso

#include "math.hpp"

#endif /* SRC_MATH2_MATH2_H_ */
