
#ifndef SRC_MATH_PRIMITIVES_MATRIX_CSC_H_
#define SRC_MATH_PRIMITIVES_MATRIX_CSC_H_

#include "matrix_csr.h"
#include "basis/containers/allocators.h"
#include "esinfo/eslog.hpp"

namespace espreso {

struct Matrix_CSC_Solver;

template <typename T, typename I>
struct _Matrix_CSC {
	I nrows, ncols, nnz, *rows, *cols;
	T *vals;
};

template <typename T, typename I = int, typename A = cpu_allocator>
class Matrix_CSC: public _Matrix_CSC<T, I>
{
public:
	Matrix_CSC(const A &ator_ = A()): _Matrix_CSC<T, I>{}, type{Matrix_Type::REAL_STRUCTURALLY_SYMMETRIC}, shape{Matrix_Shape::FULL}, _allocated{}, ator(ator_)
	{

	}

	Matrix_CSC(const Matrix_CSC &other): _Matrix_CSC<T, I>{}, type{Matrix_Type::REAL_STRUCTURALLY_SYMMETRIC}, shape{Matrix_Shape::FULL}, _allocated{}, ator(other.ator)
	{
		type = other.type;
		shape = other.shape;
		realloc(_allocated, other.nrows, other.ncols, other.nnz);
		_Matrix_CSC<T, I>::operator=(_allocated);
		for (I i = 0; i <= other.nrows; ++i) {
			this->rows[i] = other.rows[i];
		}
		for (I i = 0; i < other.nnz; ++i) {
			this->cols[i] = other.cols[i];
			this->vals[i] = other.vals[i];
		}
	}

	Matrix_CSC(Matrix_CSC &&other): _Matrix_CSC<T, I>{}, type{Matrix_Type::REAL_STRUCTURALLY_SYMMETRIC}, shape{Matrix_Shape::FULL}, _allocated{}, ator(std::move(other.ator))
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

	void resize(I nrows, I ncols, I nnz)
	{
		realloc(_allocated, nrows, ncols, nnz);
		_Matrix_CSC<T, I>::operator=(_allocated);
	}

	template<typename T2, typename I2, typename A2>
	void resize(const Matrix_CSC<T2,I2,A2> &other)
	{
		resize(other.nrows, other.ncols, other.nnz);
	}

	template<typename T2>
	void pattern(const Matrix_CSC<T2,I,A> &other)
	{
		if constexpr(!A::always_equal) if(this->ator != other.ator) eslog::error("not implemented for unequal allocators");
		realloc(_allocated, other);
		_Matrix_CSC<T, I>::operator=(_allocated);
		this->rows = other.rows;
		this->cols = other.cols;
	}

	void shallowCopy(const Matrix_CSC &other)
	{
		if constexpr(!A::always_equal) if(this->ator != other.ator) eslog::error("not implemented for unequal allocators");
		type = other.type;
		shape = other.shape;
		_Matrix_CSC<T, I>::operator=(other);
	}

	Matrix_Type type;
	Matrix_Shape shape;
	_Matrix_CSC<T, I> _allocated;
	A ator;

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

	void realloc(_Matrix_CSC<T, I> &m, I nrows, I ncols, I nnz)
	{
		if (m.nrows < nrows) {
			if (m.cols) { ator.deallocate(m.cols); m.cols = nullptr; }
			m.cols = ator.template allocate<I>(ncols + 1);
		}
		if (m.nnz < nnz) {
			if (m.rows) { ator.deallocate(m.rows); m.rows = nullptr; }
			if (m.vals) { ator.deallocate(m.vals); m.vals = nullptr; }
			m.rows = ator.template allocate<I>(nnz);
			m.vals = ator.template allocate<T>(nnz);
		}
		m.nrows = nrows;
		m.ncols = ncols;
		m.nnz = nnz;
	}

	void realloc(_Matrix_CSC<T, I> &m, const _Matrix_CSC<T, I> &other)
	{
		if (m.rows) { ator.deallocate(m.rows); m.rows = nullptr; }
		if (m.cols) { ator.deallocate(m.cols); m.cols = nullptr; }

		if (m.nnz < other.nnz) {
			clear(m);
			m.vals = ator.template allocate<T>(other.nnz);
		}
		m.nrows = other.nrows;
		m.ncols = other.ncols;
		m.nnz = other.nnz;
	}

	void clear(_Matrix_CSC<T, I> &m)
	{
		m.nrows = m.ncols = m.nnz = 0;
		if (m.rows) { ator.deallocate(m.rows); m.rows = nullptr; }
		if (m.cols) { ator.deallocate(m.cols); m.cols = nullptr; }
		if (m.vals) { ator.deallocate(m.vals); m.vals = nullptr; }
	}
};

}




#endif /* SRC_MATH_PRIMITIVES_MATRIX_CSC_H_ */
