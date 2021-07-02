
#ifndef SRC_MATH2_PRIMITIVES_MATRIX_IJV_H_
#define SRC_MATH2_PRIMITIVES_MATRIX_IJV_H_

#include "matrix_info.h"

namespace espreso {

struct Matrix_IJV_External_Representation;

struct _Matrix_IJV_Pattern {
	esint nrows, ncols, nnz, *rows, *cols;
};

template <typename T>
struct _Matrix_IJV_Vals {
	T *vals;
};

template <typename T>
struct _Matrix_IJV: public _Matrix_IJV_Pattern, public _Matrix_IJV_Vals<T> {

};

template <typename T>
class Matrix_IJV: public _Matrix_IJV<T>
{
public:
	Matrix_IJV(): _Matrix_IJV<T>{}, type{Matrix_Type::REAL_UNSYMMETRIC}, shape{Matrix_Shape::FULL}, _external{nullptr}, _allocated{}
	{

	}

	Matrix_IJV(const Matrix_IJV &other): _Matrix_IJV<T>{}, type{Matrix_Type::REAL_UNSYMMETRIC}, shape{Matrix_Shape::FULL}, _external{nullptr}, _allocated{}
	{
		this->type = other.type;
		this->shape = other.shape;
		realloc(_allocated, other.nrows, other.ncols, other.nnz);
		_Matrix_IJV<T>::operator=(_allocated);
		for (esint i = 0; i < other.nnz; ++i) {
			this->rows[i] = other.rows[i];
			this->cols[i] = other.cols[i];
			this->vals[i] = other.vals[i];
		}
	}

	Matrix_IJV(Matrix_IJV &&other): _Matrix_IJV<T>{}, type{Matrix_Type::REAL_UNSYMMETRIC}, shape{Matrix_Shape::FULL}, _external{nullptr}, _allocated{}
	{
		this->type = other.type;
		this->shape = other.shape;
		swap(*this, other);
		swap(_allocated, other._allocated);
	}

	Matrix_IJV& operator=(const Matrix_IJV &other) = delete;
//	{
//		this->type = other.type;
//		this->shape = other.shape;
//		realloc(_allocated, other.nrows, other.ncols, other.nnz);
//		_Matrix_IJV<T>::operator=(_allocated);
//		for (esint i = 0; i < other.nnz; ++i) {
//			this->rows[i] = other.rows[i];
//			this->cols[i] = other.cols[i];
//			this->vals[i] = other.vals[i];
//		}
//		return *this;
//	}

	Matrix_IJV& operator=(Matrix_IJV &&other) = delete;
//	{
//		this->type = other.type;
//		this->shape = other.shape;
//		swap(*this, other);
//		swap(_allocated, other._allocated);
//		return *this;
//	}

	~Matrix_IJV()
	{
		clear(_allocated);
	}

	void resize(esint nrows, esint ncols, esint nnz)
	{
		realloc(_allocated, nrows, ncols, nnz);
		_Matrix_IJV<T>::operator=(_allocated);
	}

	void resize(const Matrix_IJV &other)
	{
		resize(other.nrows, other.ncols, other.nnz);
	}

	void pattern(const Matrix_IJV &other)
	{
		realloc(_allocated, other);
		_Matrix_IJV<T>::operator=(_allocated);
		this->rows = other.rows;
		this->cols = other.cols;
	}

	Matrix_Type type;
	Matrix_Shape shape;
	Matrix_IJV_External_Representation *_external;

protected:
	template <typename Type>
	void swap(Type &v, Type &u)
	{
		Type tmp = v; v = u; u = tmp;
	}

	void swap(_Matrix_IJV<T> &m, _Matrix_IJV<T> &n)
	{
		swap(m.nrows, n.nrows);
		swap(m.ncols, n.ncols);
		swap(m.nnz, n.nnz);
		swap(m.rows, n.rows);
		swap(m.cols, n.cols);
		swap(m.vals, n.vals);
	}

	void realloc(_Matrix_IJV<T> &m, esint nrows, esint ncols, esint nnz)
	{
		if (m.nnz < nnz) {
			if (m.rows) { delete[] m.rows; m.rows = nullptr; }
			if (m.cols) { delete[] m.cols; m.cols = nullptr; }
			if (m.vals) { delete[] m.vals; m.vals = nullptr; }
			m.rows = new esint[nnz];
			m.cols = new esint[nnz];
			m.vals = new T[nnz];
		}
		m.nrows = nrows;
		m.ncols = ncols;
		m.nnz = nnz;
	}

	void realloc(_Matrix_IJV<T> &m, const _Matrix_IJV_Pattern &other)
	{
		if (m.rows) { delete[] m.rows; m.rows = nullptr; }
		if (m.cols) { delete[] m.cols; m.cols = nullptr; }

		if (m.nnz < other.nnz) {
			clear(m);
			m.vals = new T[other.nnz];
		}
		m.nrows = other.nrows;
		m.ncols = other.ncols;
		m.nnz = other.nnz;
	}

	void clear(_Matrix_IJV<T> &m)
	{
		m.nrows = m.ncols = m.nnz = 0;
		if (m.rows) { delete[] m.rows; m.rows = nullptr; }
		if (m.cols) { delete[] m.cols; m.cols = nullptr; }
		if (m.vals) { delete[] m.vals; m.vals = nullptr; }
	}

	_Matrix_IJV<T> _allocated;
};

}

#endif /* SRC_MATH2_PRIMITIVES_MATRIX_IJV_H_ */
