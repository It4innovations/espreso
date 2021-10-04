
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

	// y = alpha * A * x + beta * y
	template <typename T> void apply(Vector_Dense<T> &y, const T &alpha, Matrix_Dense<T> &a, const T &beta, const Vector_Dense<T> &x);
	template <typename T> void apply(Vector_Dense<T> &y, const T &alpha, Matrix_CSR<T>   &a, const T &beta, const Vector_Dense<T> &x);
	template <typename T> void apply(Vector_Dense<T> &y, const T &alpha, Matrix_IJV<T>   &a, const T &beta, const Vector_Dense<T> &x);

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

	template <typename T> void multiplyPattern(Vector_Dense<T>  &x, const Vector_Dense<T>  &y, int size, int multipicity);
	template <typename T> void multiplyPattern(Vector_Sparse<T> &x, const Vector_Sparse<T> &y, int size, int multipicity);
	template <typename T> void multiplyPattern(Matrix_Dense<T>  &x, const Matrix_Dense<T>  &y, int size, int multipicity);
	template <typename T> void multiplyPattern(Matrix_CSR<T>    &x, const Matrix_CSR<T>    &y, int size, int multipicity);
	template <typename T> void multiplyPattern(Matrix_IJV<T>    &x, const Matrix_IJV<T>    &y, int size, int multipicity);

	template <typename T> void set(Vector_Dense<T>  &x, const T &value) { set(x.size           , x.vals, 1, value); }
	template <typename T> void set(Vector_Sparse<T> &x, const T &value) { set(x.nnz            , x.vals, 1, value); }
	template <typename T> void set(Matrix_Dense<T>  &x, const T &value) { set(x.nrows * x.ncols, x.vals, 1, value); }
	template <typename T> void set(Matrix_CSR<T>    &x, const T &value) { set(x.nnz            , x.vals, 1, value); }
	template <typename T> void set(Matrix_IJV<T>    &x, const T &value) { set(x.nnz            , x.vals, 1, value); }

	// type to type
	template <typename T> void copy(Vector_Dense<T>  &x, const Vector_Dense<T>  &y) { copy(x.size           , x.vals, 1, y.vals, 1); }
	template <typename T> void copy(Vector_Sparse<T> &x, const Vector_Sparse<T> &y) { copy(x.nnz            , x.vals, 1, y.vals, 1); }
	template <typename T> void copy(Matrix_Dense<T>  &x, const Matrix_Dense<T>  &y) { copy(x.nrows * x.ncols, x.vals, 1, y.vals, 1); }
	template <typename T> void copy(Matrix_CSR<T>    &x, const Matrix_CSR<T>    &y) { copy(x.nnz            , x.vals, 1, y.vals, 1); }
	template <typename T> void copy(Matrix_IJV<T>    &x, const Matrix_IJV<T>    &y) { copy(x.nnz            , x.vals, 1, y.vals, 1); }

	// type converters
	template <typename T> void copy(Matrix_Dense<T>  &x, const Matrix_CSR<T>   &y);
	template <typename T> void copy(Matrix_Dense<T>  &x, const Matrix_IJV<T>   &y);
	template <typename T> void copy(Matrix_CSR<T>    &x, const Matrix_Dense<T> &y);
	template <typename T> void copy(Matrix_CSR<T>    &x, const Matrix_IJV<T>   &y);
	template <typename T> void copy(Matrix_IJV<T>    &x, const Matrix_Dense<T> &y);
	template <typename T> void copy(Matrix_IJV<T>    &x, const Matrix_CSR<T>   &y);

	// type to type with offsets
	template <typename T> void copy(Vector_Dense<T>  &x, const Vector_Dense<T>  &y, int offset, int size, int step);
	template <typename T> void copy(Vector_Sparse<T> &x, const Vector_Sparse<T> &y, int offset, int size, int step);
	template <typename T> void copy(Vector_Dense<T>  &x, const Vector_Sparse<T> &y, int offset, int size, int step);
	template <typename T> void copy(Vector_Sparse<T> &x, const Vector_Dense<T>  &y, int offset, int size, int step);
	template <typename T> void copy(Matrix_Dense<T>  &x, const Matrix_Dense<T>  &y, int rowOffset, int colOffset, int size, int step);
	template <typename T> void copy(Matrix_Dense<T>  &x, const Matrix_CSR<T>    &y, int rowOffset, int colOffset, int size, int step);
	template <typename T> void copy(Matrix_Dense<T>  &x, const Matrix_IJV<T>    &y, int rowOffset, int colOffset, int size, int step);
	template <typename T> void copy(Matrix_CSR<T>    &x, const Matrix_Dense<T>  &y, int rowOffset, int colOffset, int size, int step);
	template <typename T> void copy(Matrix_CSR<T>    &x, const Matrix_CSR<T>    &y, int rowOffset, int colOffset, int size, int step);
	template <typename T> void copy(Matrix_CSR<T>    &x, const Matrix_IJV<T>    &y, int rowOffset, int colOffset, int size, int step);
	template <typename T> void copy(Matrix_IJV<T>    &x, const Matrix_Dense<T>  &y, int rowOffset, int colOffset, int size, int step);
	template <typename T> void copy(Matrix_IJV<T>    &x, const Matrix_CSR<T>    &y, int rowOffset, int colOffset, int size, int step);
	template <typename T> void copy(Matrix_IJV<T>    &x, const Matrix_IJV<T>    &y, int rowOffset, int colOffset, int size, int step);

	// real to complex
	template <typename T> void copy(Matrix_CSR<std::complex<T> >    &x, const int offsetX, const Matrix_CSR<T>    &y) { copy(x.nnz, reinterpret_cast<T*>(x.vals) + offsetX, 2, y.vals, 1); }
	template <typename T> void copy(Vector_Dense<std::complex<T> >  &x, const int offsetX, const Vector_Dense<T>  &y) { copy(x.size, reinterpret_cast<T*>(x.vals) + offsetX, 2, y.vals, 1); }
	template <typename T> void copy(Vector_Sparse<std::complex<T> > &x, const int offsetX, const Vector_Sparse<T> &y) { copy(x.nnz , reinterpret_cast<T*>(x.vals) + offsetX, 2, y.vals, 1); }
	template <typename T> void copy(Vector_Dense<std::complex<T> >  &x, const int offsetX, const Vector_Sparse<T> &y) { eslog::error("call empty function getDiagonal\n"); }
	template <typename T> void copy(Vector_Sparse<std::complex<T> > &x, const int offsetX, const Vector_Dense<T>  &y) { eslog::error("call empty function getDiagonal\n"); }

	// complex to real
	template <typename T> void copy(Vector_Dense<T >  &x, const Vector_Dense<std::complex<T> >  &y, const int offsetY) { copy(x.size, x.vals, 1, reinterpret_cast<T*>(y.vals) + offsetY, 2); }
	template <typename T> void copy(Vector_Sparse<T > &x, const Vector_Sparse<std::complex<T> > &y, const int offsetY) { copy(x.nnz , x.vals, 1, reinterpret_cast<T*>(y.vals) + offsetY, 2); }
	template <typename T> void copy(Vector_Dense<T >  &x, const Vector_Sparse<std::complex<T> > &y, const int offsetY) { eslog::error("call empty function getDiagonal\n"); }
	template <typename T> void copy(Vector_Sparse<T > &x, const Vector_Dense<std::complex<T> >  &y, const int offsetY) { eslog::error("call empty function getDiagonal\n"); }

	template <typename T> void copy(Vector_Dense<T> &x, const Vector_Sparse<T> &y)
	{
		for (esint i = 0; i < y.nnz; ++i) {
			x.vals[y.indices[i]] = y.vals[i];
		}
	}
	template <typename T> void copy(Vector_Sparse<T> &x, const Vector_Dense<T> &y)
	{
		for (esint i = 0; i < x.nnz; ++i) {
			x.vals[i] = y.vals[x.indices[i]];
		}
	}

	template <typename T> void scale(const T &alpha, Vector_Dense<T>  &x) { scale(x.size           , alpha, x.vals, 1); }
	template <typename T> void scale(const T &alpha, Vector_Sparse<T> &x) { scale(x.nnz            , alpha, x.vals, 1); }
	template <typename T> void scale(const T &alpha, Matrix_Dense<T>  &x) { scale(x.nrows * x.ncols, alpha, x.vals, 1); }
	template <typename T> void scale(const T &alpha, Matrix_CSR<T>    &x) { scale(x.nnz            , alpha, x.vals, 1); }
	template <typename T> void scale(const T &alpha, Matrix_IJV<T>    &x) { scale(x.nnz            , alpha, x.vals, 1); }

	// x = alpha * y
	template <typename T> void add(Vector_Dense<T>  &x, const T &alpha, const Vector_Dense<T>  &y) { add(x.size           , x.vals, 1, alpha, y.vals, 1); }
	template <typename T> void add(Vector_Sparse<T> &x, const T &alpha, const Vector_Sparse<T> &y) { add(x.nnz            , x.vals, 1, alpha, y.vals, 1); }
	template <typename T> void add(Matrix_Dense<T>  &x, const T &alpha, const Matrix_Dense<T>  &y) { add(x.nrows * x.ncols, x.vals, 1, alpha, y.vals, 1); }
	template <typename T> void add(Matrix_CSR<T>    &x, const T &alpha, const Matrix_CSR<T>    &y) { add(x.nnz            , x.vals, 1, alpha, y.vals, 1); }
	template <typename T> void add(Matrix_IJV<T>    &x, const T &alpha, const Matrix_IJV<T>    &y) { add(x.nnz            , x.vals, 1, alpha, y.vals, 1); }
	template <typename T> void add(Vector_Sparse<T> &x, const T &alpha, const Vector_Dense<T>  &y)
	{
		for (esint i = 0; i < x.nnz; ++i) {
			x.vals[i] += alpha * y.vals[x.indices[i]];
		}
	}
	template <typename T> void add(Vector_Dense<T> &x, const T &alpha, const Vector_Sparse<T>  &y)
	{
		for (esint i = 0; i < y.nnz; ++i) {
			x.vals[y.indices[i]] += alpha * y.vals[i];
		}
	}

	template <typename T> void add(Matrix_Dense<T>  &x, const T &alpha, const Matrix_CSR<T>   &y);
	template <typename T> void add(Matrix_Dense<T>  &x, const T &alpha, const Matrix_IJV<T>   &y);
	template <typename T> void add(Matrix_CSR<T>    &x, const T &alpha, const Matrix_Dense<T> &y);
	template <typename T> void add(Matrix_CSR<T>    &x, const T &alpha, const Matrix_IJV<T>   &y);
	template <typename T> void add(Matrix_IJV<T>    &x, const T &alpha, const Matrix_Dense<T> &y);
	template <typename T> void add(Matrix_IJV<T>    &x, const T &alpha, const Matrix_CSR<T>   &y);

	// x += alpha * y [.. offset -> size  / step]
	template <typename T> void add(Vector_Sparse<std::complex<T>> &x, const int offsetX, const T &alpha, const Vector_Sparse<T> &y) { add(x.nnz, reinterpret_cast<T*>(x.vals) + offsetX, 2, alpha, y.vals, 1); }
	template <typename T> void add(Vector_Dense<std::complex<T>>  &x, const int offsetX, const T &alpha, const Vector_Dense<T>  &y) { add(x.size, reinterpret_cast<T*>(x.vals) + offsetX, 2, alpha, y.vals, 1); }
	template <typename T> void add(Matrix_CSR<std::complex<T>>    &x, const int offsetX, const T &alpha, const Matrix_CSR<T>    &y) { add(x.nnz, reinterpret_cast<T*>(x.vals) + offsetX, 2, alpha, y.vals, 1); }

	template <typename T> void add(Vector_Dense<T>  &x, const T &alpha, const Vector_Dense<T>  &y, int offset, int size, int step);
	template <typename T> void add(Vector_Sparse<T> &x, const T &alpha, const Vector_Sparse<T> &y, int offset, int size, int step);
	template <typename T> void add(Vector_Sparse<T> &x, const T &alpha, const Vector_Dense<T>  &y, int offset, int size, int step) {
		eslog::error("call empty function add\n");
	}
	template <typename T> void add(Vector_Dense<T>  &x, const T &alpha, const Vector_Sparse<T> &y, int offset, int size, int step) {
		eslog::error("call empty function add\n");
	}

	template <typename T> void add(Matrix_Dense<T>  &x, const T &alpha, const Matrix_Dense<T>  &y, int rowOffset, int colOffset, int size, int step);
	template <typename T> void add(Matrix_CSR<T>    &x, const T &alpha, const Matrix_CSR<T>    &y, int rowOffset, int colOffset, int size, int step);
	template <typename T> void add(Matrix_IJV<T>    &x, const T &alpha, const Matrix_IJV<T>    &y, int rowOffset, int colOffset, int size, int step);

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

	template <typename T> void getDiagonal(const Matrix_CSR<T> &m, Vector_Dense<T> &v)
	{
		eslog::error("call empty function getDiagonal\n");
	}

	template <class T> void store(const T &x, const char* file);

} // math
} // espreso

#include "math2.hpp"

#endif /* SRC_MATH2_MATH2_H_ */
