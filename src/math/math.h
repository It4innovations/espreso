
#ifndef SRC_MATH_MATH_H_
#define SRC_MATH_MATH_H_

#include "esinfo/eslog.h"
#include "primitives/vector_dense.h"
#include "primitives/vector_sparse.h"
#include "primitives/matrix_dense.h"
#include "primitives/matrix_csr.h"
#include "primitives/matrix_ijv.h"

#include "wrappers/math.blas.h"
#include "wrappers/math.spblas.h"
#include "wrappers/math.lapack.h"
#include "wrappers/math.spsolver.h"

#include <complex>

namespace espreso {
namespace math { // interface to wrappers

	// x = value
	template <typename T, typename I>
	void set(const I size, T *x, const int incX, const T &value)
	{
		for (esint i = 0; i < size; i += incX) {
			x[i] = value;
		}
	}

} // math (interface to wrappers)
} // espreso

namespace espreso {
namespace math {

	template <typename T, typename I> void copy(Vector_Dense<T, I>  &x, const Vector_Dense<T, I>  &y) { blas::copy(x.size           , x.vals, 1, y.vals, 1); }
	template <typename T, typename I> void copy(Vector_Sparse<T, I> &x, const Vector_Sparse<T, I> &y) { blas::copy(x.nnz            , x.vals, 1, y.vals, 1); }
	template <typename T, typename I> void copy(Matrix_Dense<T, I>  &x, const Matrix_Dense<T, I>  &y) { blas::copy(x.nrows * x.ncols, x.vals, 1, y.vals, 1); }
	template <typename T, typename I> void copy(Matrix_CSR<T, I>    &x, const Matrix_CSR<T, I>    &y) { blas::copy(x.nnz            , x.vals, 1, y.vals, 1); }
	template <typename T, typename I> void copy(Matrix_IJV<T, I>    &x, const Matrix_IJV<T, I>    &y) { blas::copy(x.nnz            , x.vals, 1, y.vals, 1); }

	template <typename T, typename I> void scale(const T &alpha, Vector_Dense<T, I>  &x) { blas::scale(x.size           , alpha, x.vals, 1); }
	template <typename T, typename I> void scale(const T &alpha, Vector_Sparse<T, I> &x) { blas::scale(x.nnz            , alpha, x.vals, 1); }
	template <typename T, typename I> void scale(const T &alpha, Matrix_Dense<T, I>  &x) { blas::scale(x.nrows * x.ncols, alpha, x.vals, 1); }
	template <typename T, typename I> void scale(const T &alpha, Matrix_CSR<T, I>    &x) { blas::scale(x.nnz            , alpha, x.vals, 1); }
	template <typename T, typename I> void scale(const T &alpha, Matrix_IJV<T, I>    &x) { blas::scale(x.nnz            , alpha, x.vals, 1); }

	template <typename T, typename I> void add(Vector_Dense<T, I>  &x, const T &alpha, const Vector_Dense<T, I>  &y) { blas::add(x.size           , x.vals, 1, alpha, y.vals, 1); }
	template <typename T, typename I> void add(Vector_Sparse<T, I> &x, const T &alpha, const Vector_Sparse<T, I> &y) { blas::add(x.nnz            , x.vals, 1, alpha, y.vals, 1); }
	template <typename T, typename I> void add(Matrix_Dense<T, I>  &x, const T &alpha, const Matrix_Dense<T, I>  &y) { blas::add(x.nrows * x.ncols, x.vals, 1, alpha, y.vals, 1); }
	template <typename T, typename I> void add(Matrix_CSR<T, I>    &x, const T &alpha, const Matrix_CSR<T, I>    &y) { blas::add(x.nnz            , x.vals, 1, alpha, y.vals, 1); }
	template <typename T, typename I> void add(Matrix_IJV<T, I>    &x, const T &alpha, const Matrix_IJV<T, I>    &y) { blas::add(x.nnz            , x.vals, 1, alpha, y.vals, 1); }

	template <typename T, typename I> T dot(const Vector_Dense<T, I>  &x, const Vector_Dense<T, I>  &y) { return blas::dot(x.size, x.vals, 1, y.vals, 1); }
	template <typename T, typename I> T dot(const Vector_Sparse<T, I> &x, const Vector_Sparse<T, I> &y) { return blas::dot(x.nnz , x.vals, 1, y.vals, 1); }

	template <typename T, typename I> T norm(const Vector_Dense<T, I>  &x) { return blas::norm(x.size, x.vals, 1); }
	template <typename T, typename I> T norm(const Vector_Sparse<T, I> &x) { return blas::norm(x.nnz , x.vals, 1); }

	template <typename T, typename I> void set(Vector_Dense<T, I>  &x, const T &value) { set(x.size           , x.vals, 1, value); }
	template <typename T, typename I> void set(Vector_Sparse<T, I> &x, const T &value) { set(x.nnz            , x.vals, 1, value); }
	template <typename T, typename I> void set(Matrix_Dense<T, I>  &x, const T &value) { set(x.nrows * x.ncols, x.vals, 1, value); }
	template <typename T, typename I> void set(Matrix_CSR<T, I>    &x, const T &value) { set(x.nnz            , x.vals, 1, value); }
	template <typename T, typename I> void set(Matrix_IJV<T, I>    &x, const T &value) { set(x.nnz            , x.vals, 1, value); }

	template <typename T, typename I> void combine(Matrix_CSR<T, I> &C, const Matrix_CSR<T, I> &A, const Matrix_CSR<T, I> &B)
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

	template <typename T, typename I> void sumCombined(Matrix_CSR<T, I> &C, const T &alpha, const Matrix_CSR<T, I> &A, const Matrix_CSR<T, I> &B)
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
			T      *a = A.vals + A.rows[r    ] - Indexing::CSR;
			T      *b = B.vals + B.rows[r    ] - Indexing::CSR;
			T      *c = C.vals + C.rows[r    ] - Indexing::CSR;
			while (bC != eC) {
				*c = 0;
				if (bA != eA && *bC == *bA) { *c += *a++; ++bA; }
				if (bB != eB && *bC == *bB) { *c += *b++; ++bB; }
				++bC; ++c;
			}
		}
	}

	template <typename T, typename I> void orthonormalize(Matrix_Dense<T, I> &m)
	{
		for (esint r = 0; r < m.nrows; ++r) {
			for (esint rr = 0; rr < r; ++rr) {
				T scale = math::blas::dot(m.ncols, m.vals + rr * m.ncols, 1, m.vals + r * m.ncols, 1) / math::blas::dot(m.ncols, m.vals + rr * m.ncols, 1, m.vals + rr * m.ncols, 1);
				math::blas::add(m.ncols, m.vals + r * m.ncols, 1, -scale, m.vals + rr * m.ncols, 1);
			}
			math::blas::scale(m.ncols, T{1.} / math::blas::norm(m.ncols, m.vals + r * m.ncols, 1), m.vals + r * m.ncols, 1);
		}
	}

	template <class T> void store(const T &x, const char* file);

} // math
} // espreso

#endif /* SRC_MATH_MATH_H_ */
