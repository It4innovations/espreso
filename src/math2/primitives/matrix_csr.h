
#ifndef SRC_MATH2_PRIMITIVES_MATRIX_CSR_H_
#define SRC_MATH2_PRIMITIVES_MATRIX_CSR_H_

#include "matrix_info.h"

namespace espreso {

struct Matrix_CSR_External_Representation;

struct _Matrix_CSR_Pattern {
	enum: int {
		Indexing = 1
	};

	esint nrows, ncols, nnz, *rows, *cols;
};

template <typename T>
struct _Matrix_CSR_Vals {
	T *vals;
};

template <typename T>
struct _Matrix_CSR: public _Matrix_CSR_Pattern, public _Matrix_CSR_Vals<T> {

};

template <typename T>
class Matrix_CSR: public _Matrix_CSR<T>
{
public:
	Matrix_CSR(): _Matrix_CSR<T>{}, type{Matrix_Type::REAL_UNSYMMETRIC}, shape{Matrix_Shape::FULL}, _external{nullptr}, _allocated{}
	{

	}

	Matrix_CSR(const Matrix_CSR &other): _Matrix_CSR<T>{}, type{Matrix_Type::REAL_UNSYMMETRIC}, shape{Matrix_Shape::FULL}, _external{nullptr}, _allocated{}
	{
		type = other.type;
		shape = other.shape;
		realloc(_allocated, other.nrows, other.ncols, other.nnz);
		_Matrix_CSR<T>::operator=(_allocated);
		for (esint i = 0; i <= other.nrows; ++i) {
			this->rows[i] = other.rows[i];
		}
		for (esint i = 0; i < other.nnz; ++i) {
			this->cols[i] = other.cols[i];
			this->vals[i] = other.vals[i];
		}
	}

	Matrix_CSR(Matrix_CSR &&other): _Matrix_CSR<T>{}, type{Matrix_Type::REAL_UNSYMMETRIC}, shape{Matrix_Shape::FULL}, _external{nullptr}, _allocated{}
	{
		type = other.type;
		shape = other.shape;
		swap(_external, other._external);
		swap(*this, other);
		swap(_allocated, other._allocated);
	}

	Matrix_CSR& operator=(const Matrix_CSR &other) = delete;
//	{
//		type = other.type;
//		shape = other.shape;
//		realloc(_allocated, other.nrows, other.ncols, other.nnz);
//		_Matrix_CSR<T>::operator=(_allocated);
//		for (esint i = 0; i <= other.nrows; ++i) {
//			this->rows[i] = other.rows[i];
//		}
//		for (esint i = 0; i < other.nnz; ++i) {
//			this->cols[i] = other.cols[i];
//			this->vals[i] = other.vals[i];
//		}
//		return *this;
//	}

	Matrix_CSR& operator=(Matrix_CSR &&other) = delete;
//	{
//		type = other.type;
//		shape = other.shape;
//		swap(_external, other._external);
//		swap(*this, other);
//		swap(_allocated, other._allocated);
//		return *this;
//	}

	~Matrix_CSR()
	{
		clear(_allocated);
	}

	void resize(esint nrows, esint ncols, esint nnz)
	{
		realloc(_allocated, nrows, ncols, nnz);
		_Matrix_CSR<T>::operator=(_allocated);
	}

	void resize(const Matrix_CSR &other)
	{
		resize(other.nrows, other.ncols, other.nnz);
	}

	void pattern(const Matrix_CSR &other)
	{
		realloc(_allocated, other);
		_Matrix_CSR<T>::operator=(_allocated);
		this->rows = other.rows;
		this->cols = other.cols;
	}

	Matrix_Type type;
	Matrix_Shape shape;
	Matrix_CSR_External_Representation *_external;

protected:
	template <typename Type>
	void swap(Type &v, Type &u)
	{
		Type tmp = v; v = u; u = tmp;
	}

	void swap(_Matrix_CSR<T> &m, _Matrix_CSR<T> &n)
	{
		swap(m.nrows, n.nrows);
		swap(m.ncols, n.ncols);
		swap(m.nnz, n.nnz);
		swap(m.rows, n.rows);
		swap(m.cols, n.cols);
		swap(m.vals, n.vals);
	}

	void realloc(_Matrix_CSR<T> &m, esint nrows, esint ncols, esint nnz)
	{
		if (m.nrows < nrows) {
			if (m.rows) { delete[] m.rows; m.rows = nullptr; }
			m.rows = new esint[nrows + 1];
		}
		if (m.nnz < nnz) {
			if (m.cols) { delete[] m.cols; m.cols = nullptr; }
			if (m.vals) { delete[] m.vals; m.vals = nullptr; }
			m.cols = new esint[nnz];
			m.vals = new T[nnz];
		}
		m.nrows = nrows;
		m.ncols = ncols;
		m.nnz = nnz;
	}

	void realloc(_Matrix_CSR<T> &m, const _Matrix_CSR_Pattern &other)
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

	void clear(_Matrix_CSR<T> &m)
	{
		m.nrows = m.ncols = m.nnz = 0;
		if (m.rows) { delete[] m.rows; m.rows = nullptr; }
		if (m.cols) { delete[] m.cols; m.cols = nullptr; }
		if (m.vals) { delete[] m.vals; m.vals = nullptr; }
	}

	_Matrix_CSR<T> _allocated;
};

}

#endif /* SRC_MATH2_PRIMITIVES_MATRIX_CSR_H_ */
