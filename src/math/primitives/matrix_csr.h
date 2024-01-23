
#ifndef SRC_MATH2_PRIMITIVES_MATRIX_CSR_H_
#define SRC_MATH2_PRIMITIVES_MATRIX_CSR_H_

#include "matrix_info.h"
#include "basis/containers/allocators.h"

namespace espreso {

struct Matrix_CSR_Solver;

template <typename T, typename I>
struct _Matrix_CSR {
	I nrows, ncols, nnz, *rows, *cols;
	T *vals;
};

template <typename T, typename I = int, template<typename> typename A = cpu_allocator>
class Matrix_CSR: public _Matrix_CSR<T, I>
{
public:
	Matrix_CSR(): _Matrix_CSR<T, I>{}, type{Matrix_Type::REAL_STRUCTURALLY_SYMMETRIC}, shape{Matrix_Shape::FULL}, _allocated{}
	{

	}

	Matrix_CSR(const Matrix_CSR &other): _Matrix_CSR<T, I>{}, type{Matrix_Type::REAL_STRUCTURALLY_SYMMETRIC}, shape{Matrix_Shape::FULL}, _allocated{}
	{
		type = other.type;
		shape = other.shape;
		realloc(_allocated, other.nrows, other.ncols, other.nnz);
		_Matrix_CSR<T, I>::operator=(_allocated);
		for (esint i = 0; i <= other.nrows; ++i) {
			this->rows[i] = other.rows[i];
		}
		for (esint i = 0; i < other.nnz; ++i) {
			this->cols[i] = other.cols[i];
			this->vals[i] = other.vals[i];
		}
	}

	Matrix_CSR(Matrix_CSR &&other): _Matrix_CSR<T, I>{}, type{Matrix_Type::REAL_STRUCTURALLY_SYMMETRIC}, shape{Matrix_Shape::FULL}, _allocated{}
	{
		swap(*static_cast<_Matrix_CSR<T, I>*>(this), *static_cast<_Matrix_CSR<T, I>*>(&other));
		type = other.type;
		shape = other.shape;
		swap(_allocated, other._allocated);
	}

	Matrix_CSR& operator=(const Matrix_CSR &other) = delete;
//	{
//		type = other.type;
//		shape = other.shape;
//		realloc(_allocated, other.nrows, other.ncols, other.nnz);
//		_Matrix_CSR<T, I>::operator=(_allocated);
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
		_Matrix_CSR<T, I>::operator=(_allocated);
	}

	void resize(const Matrix_CSR &other)
	{
		resize(other.nrows, other.ncols, other.nnz);
	}

	void pattern(const Matrix_CSR &other)
	{
		realloc(_allocated, other);
		_Matrix_CSR<T, I>::operator=(_allocated);
		this->rows = other.rows;
		this->cols = other.cols;
	}

	void shallowCopy(const Matrix_CSR &other)
	{
		type = other.type;
		shape = other.shape;
		_Matrix_CSR<T, I>::operator=(other);
	}

	Matrix_Type type;
	Matrix_Shape shape;
	_Matrix_CSR<T, I> _allocated;

protected:
	template <typename Type>
	void swap(Type &v, Type &u)
	{
		Type tmp = v; v = u; u = tmp;
	}

	void swap(_Matrix_CSR<T, I> &m, _Matrix_CSR<T, I> &n)
	{
		swap(m.nrows, n.nrows);
		swap(m.ncols, n.ncols);
		swap(m.nnz, n.nnz);
		swap(m.rows, n.rows);
		swap(m.cols, n.cols);
		swap(m.vals, n.vals);
	}

	void realloc(_Matrix_CSR<T, I> &m, esint nrows, esint ncols, esint nnz)
	{
		if (m.nrows < nrows) {
			if (m.rows) { delete[] m.rows; m.rows = nullptr; }
		}
		if (m.rows == nullptr) {
			m.rows = new esint[nrows + 1];
		}
		if (m.nnz < nnz) {
			if (m.cols) { delete[] m.cols; m.cols = nullptr; }
			if (m.vals) { delete[] m.vals; m.vals = nullptr; }
		}
		if (m.cols == nullptr) {
			m.cols = new esint[nnz];
		}
		if (m.vals == nullptr) {
			m.vals = new T[nnz];
		}
		m.nrows = nrows;
		m.ncols = ncols;
		m.nnz = nnz;
	}

	void realloc(_Matrix_CSR<T, I> &m, const _Matrix_CSR<T, I> &other)
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

	void clear(_Matrix_CSR<T, I> &m)
	{
		m.nrows = m.ncols = m.nnz = 0;
		if (m.rows) { delete[] m.rows; m.rows = nullptr; }
		if (m.cols) { delete[] m.cols; m.cols = nullptr; }
		if (m.vals) { delete[] m.vals; m.vals = nullptr; }
	}
};

}

#endif /* SRC_MATH2_PRIMITIVES_MATRIX_CSR_H_ */
