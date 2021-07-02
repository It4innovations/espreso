
#ifndef SRC_MATH2_PRIMITIVES_MATRIX_DENSE_H_
#define SRC_MATH2_PRIMITIVES_MATRIX_DENSE_H_

#include "matrix_info.h"

namespace espreso {

struct Matrix_Dense_External_Representation;

struct _Matrix_Dense_Pattern {
	esint nrows, ncols;
};

template <typename T>
struct _Matrix_Dense_Vals {
	T *vals;
};

template <typename T>
struct _Matrix_Dense: public _Matrix_Dense_Pattern, public _Matrix_Dense_Vals<T> {

};

template <typename T>
class Matrix_Dense: public _Matrix_Dense<T>
{
public:
	Matrix_Dense(): _Matrix_Dense<T>{}, type{Matrix_Type::REAL_UNSYMMETRIC}, shape{Matrix_Shape::FULL}, _external{nullptr}, _allocated{}
	{

	}

	Matrix_Dense(const Matrix_Dense &other): _Matrix_Dense<T>{}, type{Matrix_Type::REAL_UNSYMMETRIC}, shape{Matrix_Shape::FULL}, _external{nullptr}, _allocated{}
	{
		this->type = other.type;
		this->shape = other.shape;
		realloc(_allocated, other.nrows, other.ncols);
		_Matrix_Dense<T>::operator=(_allocated);
		for (esint i = 0; i < other.nrows * other.ncols; ++i) {
			this->vals[i] = other.vals[i];
		}
	}

	Matrix_Dense(Matrix_Dense &&other): _Matrix_Dense<T>{}, type{Matrix_Type::REAL_UNSYMMETRIC}, shape{Matrix_Shape::FULL}, _external{nullptr}, _allocated{}
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

	void resize(esint nrows, esint ncols)
	{
		realloc(_allocated, nrows, ncols);
		_Matrix_Dense<T>::operator=(_allocated);
	}

	void resize(const Matrix_Dense &other)
	{
		resize(other.nrows, other.ncols);
	}

	void pattern(const Matrix_Dense &other)
	{
		realloc(_allocated, other.nrows, other.ncols);
		_Matrix_Dense<T>::operator=(_allocated);
	}

	Matrix_Type type;
	Matrix_Shape shape;
	Matrix_Dense_External_Representation *_external;

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

	void realloc(_Matrix_Dense<T> &m, esint nrows, esint ncols)
	{
		if (m.nrows * m.ncols < nrows * ncols) {
			clear(m);
			m.vals = new T[nrows * ncols];
		}
		m.nrows = nrows;
		m.ncols = ncols;
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
