
#ifndef SRC_MATH2_PRIMITIVES_MATRIX_DENSE_H_
#define SRC_MATH2_PRIMITIVES_MATRIX_DENSE_H_

#include "matrix_info.h"
#include "slice.h"

namespace espreso {

template <typename T>
struct _Matrix_Dense {
	esint nrows, ncols, nnz;
	T *vals;
};

template <typename T>
class Matrix_Dense: public _Matrix_Dense<T>
{
public:
	Matrix_Dense(): _Matrix_Dense<T>{}, type{Matrix_Type::REAL_STRUCTURALLY_SYMMETRIC}, shape{Matrix_Shape::FULL}, _allocated{}
	{

	}

	Matrix_Dense(const Matrix_Dense &other): _Matrix_Dense<T>{}, _allocated{}
	{
		this->type = other.type;
		this->shape = other.shape;
		realloc(_allocated, this->shape, other.nrows, other.ncols);
		_Matrix_Dense<T>::operator=(_allocated);
		for (esint i = 0; i < other.nnz; ++i) {
			this->vals[i] = other.vals[i];
		}
	}

	Matrix_Dense(Matrix_Dense &&other): _Matrix_Dense<T>{}, _allocated{}
	{
		this->type = other.type;
		this->shape = other.shape;
		swap(*this, other);
		swap(_allocated, other._allocated);
	}

	Matrix_Dense& operator=(const Matrix_Dense &other) = delete;
//	{
//		this->type = other.type;
//		this->shape = other.shape;
//		realloc(_allocated, other.nrows, other.ncols);
//		_Matrix_Dense<T>::operator=(_allocated);
//		for (esint i = 0; i < other.nrows * other.ncols; ++i) {
//			this->vals[i] = other.vals[i];
//		}
//		return *this;
//	}

	Matrix_Dense& operator=(Matrix_Dense &&other) = delete;
//	{
//		this->type = other.type;
//		this->shape = other.shape;
//		swap(*this, other);
//		swap(_allocated, other._allocated);
//		return *this;
//	}

	~Matrix_Dense()
	{
		clear(_allocated);
	}

	_Matrix_Dense<T>& allocated()
	{
		return _allocated;
	}

	void resize(esint nrows, esint ncols)
	{
		realloc(_allocated, shape, nrows, ncols);
		_Matrix_Dense<T>::operator=(_allocated);
	}

	void resize(const Matrix_Dense &other)
	{
		resize(other.nrows, other.ncols);
	}

	void pattern(const Matrix_Dense &other)
	{
		realloc(_allocated, other.shape, other.nrows, other.ncols);
		_Matrix_Dense<T>::operator=(_allocated);
	}

	void slice(const Slice &rows, const Slice &cols)
	{
		submatrix[0] = rows;
		submatrix[1] = cols;
	}

	Matrix_Type type;
	Matrix_Shape shape;
	Slice submatrix[2];

protected:
	template <typename Type>
	void swap(Type &v, Type &u)
	{
		Type tmp = v; v = u; u = tmp;
	}

	void swap(_Matrix_Dense<T> &m, _Matrix_Dense<T> &n)
	{
		swap(m.nrows, n.nrows);
		swap(m.ncols, n.ncols);
		swap(m.vals, n.vals);
	}

	void realloc(_Matrix_Dense<T> &m, const Matrix_Shape &shape, esint nrows, esint ncols)
	{
		esint nnz;
		if (shape == Matrix_Shape::FULL) {
			nnz = nrows * ncols;
		} else {
			nnz = (nrows * ncols - nrows) / 2 + nrows;
		}
		if (m.nrows * m.ncols < nrows * ncols) {
			clear(m);
			m.vals = new T[nnz];
		}
		m.nrows = nrows;
		m.ncols = ncols;
		m.nnz = nnz;
	}

	void clear(_Matrix_Dense<T> &m)
	{
		m.nrows = m.ncols = 0;
		if (m.vals) { delete[] m.vals; m.vals = nullptr; }
	}

	_Matrix_Dense<T> _allocated;
};


}

#endif /* SRC_MATH2_PRIMITIVES_MATRIX_DENSE_H_ */
