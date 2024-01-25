
#ifndef SRC_ANALYSIS_MATH_MATH_PHYSICS_H_
#define SRC_ANALYSIS_MATH_MATH_PHYSICS_H_

#include "math/math.h"

namespace espreso {
namespace math {

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
template <typename T> void copy(Matrix_Dense<std::complex<T> > &x, const int offsetX, const Matrix_CSR<T>   &y);
template <typename T> void copy(Matrix_Dense<std::complex<T> > &x, const int offsetX, const Matrix_IJV<T>   &y);
template <typename T> void copy(Matrix_CSR<std::complex<T> >   &x, const int offsetX, const Matrix_Dense<T> &y);
template <typename T> void copy(Matrix_CSR<std::complex<T> >   &x, const int offsetX, const Matrix_IJV<T>   &y);
template <typename T> void copy(Matrix_IJV<std::complex<T> >   &x, const int offsetX, const Matrix_Dense<T> &y);
template <typename T> void copy(Matrix_IJV<std::complex<T> >   &x, const int offsetX, const Matrix_CSR<T>   &y);

template <typename T> void copy(Matrix_Dense<std::complex<T> >  &x, const int offsetX, const Matrix_Dense<T>  &y) { blas::copy(x.nrows * x.ncols, reinterpret_cast<T*>(x.vals) + offsetX, 2, y.vals, 1); }
template <typename T> void copy(Matrix_CSR<std::complex<T> >    &x, const int offsetX, const Matrix_CSR<T>    &y) { blas::copy(x.nnz            , reinterpret_cast<T*>(x.vals) + offsetX, 2, y.vals, 1); }
template <typename T> void copy(Matrix_IJV<std::complex<T> >    &x, const int offsetX, const Matrix_IJV<T>    &y) { blas::copy(x.nnz            , reinterpret_cast<T*>(x.vals) + offsetX, 2, y.vals, 1); }
template <typename T> void copy(Vector_Dense<std::complex<T> >  &x, const int offsetX, const Vector_Dense<T>  &y) { blas::copy(x.size           , reinterpret_cast<T*>(x.vals) + offsetX, 2, y.vals, 1); }
template <typename T> void copy(Vector_Sparse<std::complex<T> > &x, const int offsetX, const Vector_Sparse<T> &y) { blas::copy(x.nnz            , reinterpret_cast<T*>(x.vals) + offsetX, 2, y.vals, 1); }
template <typename T> void copy(Vector_Dense<std::complex<T> >  &x, const int offsetX, const Vector_Sparse<T> &y) { eslog::error("call empty function copy\n"); }
template <typename T> void copy(Vector_Sparse<std::complex<T> > &x, const int offsetX, const Vector_Dense<T>  &y) { eslog::error("call empty function copy\n"); }

// complex to real
template <typename T> void copy(Vector_Dense<T >  &x, const Vector_Dense<std::complex<T> >  &y, const int offsetY) { blas::copy(x.size, x.vals, 1, reinterpret_cast<T*>(y.vals) + offsetY, 2); }
template <typename T> void copy(Vector_Sparse<T > &x, const Vector_Sparse<std::complex<T> > &y, const int offsetY) { blas::copy(x.nnz , x.vals, 1, reinterpret_cast<T*>(y.vals) + offsetY, 2); }
template <typename T> void copy(Vector_Dense<T >  &x, const Vector_Sparse<std::complex<T> > &y, const int offsetY) { eslog::error("call empty function copy\n"); }
template <typename T> void copy(Vector_Sparse<T > &x, const Vector_Dense<std::complex<T> >  &y, const int offsetY) { eslog::error("call empty function copy\n"); }

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


// x += alpha * y [.. offset -> size  / step]
template <typename T> void add(Matrix_Dense<std::complex<T>> &x, const int offsetX, const T &alpha, const Matrix_CSR<T>   &y);
template <typename T> void add(Matrix_Dense<std::complex<T>> &x, const int offsetX, const T &alpha, const Matrix_IJV<T>   &y);
template <typename T> void add(Matrix_CSR<std::complex<T>>   &x, const int offsetX, const T &alpha, const Matrix_Dense<T> &y);
template <typename T> void add(Matrix_CSR<std::complex<T>>   &x, const int offsetX, const T &alpha, const Matrix_IJV<T>   &y);
template <typename T> void add(Matrix_IJV<std::complex<T>>   &x, const int offsetX, const T &alpha, const Matrix_Dense<T> &y);
template <typename T> void add(Matrix_IJV<std::complex<T>>   &x, const int offsetX, const T &alpha, const Matrix_CSR<T>   &y);

template <typename T> void add(Matrix_Dense<std::complex<T>>  &x, const int offsetX, const T &alpha, const Matrix_Dense<T>  &y) { blas::add(x.nrows * x.ncols, reinterpret_cast<T*>(x.vals) + offsetX, 2, alpha, y.vals, 1); }
template <typename T> void add(Matrix_CSR<std::complex<T>>    &x, const int offsetX, const T &alpha, const Matrix_CSR<T>    &y) { blas::add(x.nnz            , reinterpret_cast<T*>(x.vals) + offsetX, 2, alpha, y.vals, 1); }
template <typename T> void add(Matrix_IJV<std::complex<T>>    &x, const int offsetX, const T &alpha, const Matrix_IJV<T>    &y) { blas::add(x.nnz            , reinterpret_cast<T*>(x.vals) + offsetX, 2, alpha, y.vals, 1); }
template <typename T> void add(Vector_Sparse<std::complex<T>> &x, const int offsetX, const T &alpha, const Vector_Sparse<T> &y) { blas::add(x.nnz            , reinterpret_cast<T*>(x.vals) + offsetX, 2, alpha, y.vals, 1); }
template <typename T> void add(Vector_Dense<std::complex<T>>  &x, const int offsetX, const T &alpha, const Vector_Dense<T>  &y) { blas::add(x.size           , reinterpret_cast<T*>(x.vals) + offsetX, 2, alpha, y.vals, 1); }

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

}
}

#include "analysis/math/math.physics.copy.hpp"
#include "analysis/math/math.physics.add.hpp"

#endif /* SRC_ANALYSIS_MATH_MATH_PHYSICS_H_ */
