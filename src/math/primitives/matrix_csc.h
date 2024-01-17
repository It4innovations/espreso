
#ifndef SRC_MATH_PRIMITIVES_MATRIX_CSC_H_
#define SRC_MATH_PRIMITIVES_MATRIX_CSC_H_

#include "matrix_csr.h"
#include "basis/containers/allocators.h"

namespace espreso {

struct Matrix_CSC_Solver;

template <typename T, typename I>
struct _Matrix_CSC {
	I nrows, ncols, nnz, *rows, *cols;
	T *vals;
};

template <typename T, typename I = int, template<typename> typename A = cpu_allocator>
class Matrix_CSC: public _Matrix_CSC<T, I>
{
public:
	Matrix_CSC(): _Matrix_CSC<T, I>{}, type{Matrix_Type::REAL_STRUCTURALLY_SYMMETRIC}, shape{Matrix_Shape::FULL}, _allocated{}
	{

	}

	Matrix_CSC(const Matrix_CSC &other): _Matrix_CSC<T, I>{}, type{Matrix_Type::REAL_STRUCTURALLY_SYMMETRIC}, shape{Matrix_Shape::FULL}, _allocated{}
	{
		type = other.type;
		shape = other.shape;
		realloc(_allocated, other.nrows, other.ncols, other.nnz);
		_Matrix_CSC<T, I>::operator=(_allocated);
		for (esint i = 0; i <= other.nrows; ++i) {
			this->rows[i] = other.rows[i];
		}
		for (esint i = 0; i < other.nnz; ++i) {
			this->cols[i] = other.cols[i];
			this->vals[i] = other.vals[i];
		}
	}

	Matrix_CSC(Matrix_CSC &&other): _Matrix_CSC<T, I>{}, type{Matrix_Type::REAL_STRUCTURALLY_SYMMETRIC}, shape{Matrix_Shape::FULL}, _allocated{}
	{
		type = other.type;
		shape = other.shape;
		swap(*this, other);
		swap(_allocated, other._allocated);
	}

	Matrix_CSC& operator=(const Matrix_CSC &other) = delete;
//	{
//		type = other.type;
//		shape = other.shape;
//		realloc(_allocated, other.nrows, other.ncols, other.nnz);
//		_Matrix_CSC<T, I>::operator=(_allocated);
//		for (esint i = 0; i <= other.nrows; ++i) {
//			this->rows[i] = other.rows[i];
//		}
//		for (esint i = 0; i < other.nnz; ++i) {
//			this->cols[i] = other.cols[i];
//			this->vals[i] = other.vals[i];
//		}
//		return *this;
//	}

	Matrix_CSC& operator=(Matrix_CSC &&other) = delete;
//	{
//		type = other.type;
//		shape = other.shape;
//		swap(_external, other._external);
//		swap(*this, other);
//		swap(_allocated, other._allocated);
//		return *this;
//	}

	~Matrix_CSC()
	{
		clear(_allocated);
	}

	void resize(esint nrows, esint ncols, esint nnz)
	{
		realloc(_allocated, nrows, ncols, nnz);
		_Matrix_CSC<T, I>::operator=(_allocated);
	}

	void resize(const Matrix_CSC &other)
	{
		resize(other.nrows, other.ncols, other.nnz);
	}

	void pattern(const Matrix_CSC &other)
	{
		realloc(_allocated, other);
		_Matrix_CSC<T, I>::operator=(_allocated);
		this->rows = other.rows;
		this->cols = other.cols;
	}

	void shallowCopy(const Matrix_CSC &other)
	{
		type = other.type;
		shape = other.shape;
		_Matrix_CSC<T, I>::operator=(other);
	}

	Matrix_Type type;
	Matrix_Shape shape;
	_Matrix_CSC<T, I> _allocated;

protected:
	template <typename Type>
	void swap(Type &v, Type &u)
	{
		Type tmp = v; v = u; u = tmp;
	}

	void swap(_Matrix_CSC<T, I> &m, _Matrix_CSC<T, I> &n)
	{
		swap(m.nrows, n.nrows);
		swap(m.ncols, n.ncols);
		swap(m.nnz, n.nnz);
		swap(m.rows, n.rows);
		swap(m.cols, n.cols);
		swap(m.vals, n.vals);
	}

	void realloc(_Matrix_CSC<T, I> &m, esint nrows, esint ncols, esint nnz)
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

	void realloc(_Matrix_CSC<T, I> &m, const _Matrix_CSC<T, I> &other)
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

	void clear(_Matrix_CSC<T, I> &m)
	{
		m.nrows = m.ncols = m.nnz = 0;
		if (m.rows) { delete[] m.rows; m.rows = nullptr; }
		if (m.cols) { delete[] m.cols; m.cols = nullptr; }
		if (m.vals) { delete[] m.vals; m.vals = nullptr; }
	}
};

}




#endif /* SRC_MATH_PRIMITIVES_MATRIX_CSC_H_ */
