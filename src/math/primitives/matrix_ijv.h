
#ifndef SRC_MATH2_PRIMITIVES_MATRIX_IJV_H_
#define SRC_MATH2_PRIMITIVES_MATRIX_IJV_H_

#include "matrix_info.h"
#include "basis/containers/allocators.h"
#include "esinfo/eslog.hpp"

namespace espreso {

template <typename T, typename I>
struct _Matrix_IJV {
	I nrows, ncols, nnz, *rows, *cols;
	T *vals;
};

template <typename T, typename I = int, typename A = cpu_allocator>
class Matrix_IJV: public _Matrix_IJV<T, I>
{
public:
	Matrix_IJV(const A &ator_ = A()): _Matrix_IJV<T, I>{}, type{Matrix_Type::REAL_NONSYMMETRIC}, shape{Matrix_Shape::FULL}, _allocated{}, ator(ator_)
	{

	}

	Matrix_IJV(const Matrix_IJV &other): _Matrix_IJV<T, I>{}, type{Matrix_Type::REAL_NONSYMMETRIC}, shape{Matrix_Shape::FULL}, _allocated{}, ator(other.ator)
	{
		this->type = other.type;
		this->shape = other.shape;
		realloc(_allocated, other.nrows, other.ncols, other.nnz);
		_Matrix_IJV<T, I>::operator=(_allocated);
		for (I i = 0; i < other.nnz; ++i) {
			this->rows[i] = other.rows[i];
			this->cols[i] = other.cols[i];
			this->vals[i] = other.vals[i];
		}
	}

	Matrix_IJV(Matrix_IJV &&other): _Matrix_IJV<T, I>{}, type{Matrix_Type::REAL_NONSYMMETRIC}, shape{Matrix_Shape::FULL}, _allocated{}, ator(std::move(other.ator))
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
//		_Matrix_IJV<T, I>::operator=(_allocated);
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

	void resize(I nrows, I ncols, I nnz)
	{
		realloc(_allocated, nrows, ncols, nnz);
		_Matrix_IJV<T, I>::operator=(_allocated);
	}

	template<typename T2, typename I2, typename A2>
	void resize(const Matrix_IJV<T2,I2,A2> &other)
	{
		resize(other.nrows, other.ncols, other.nnz);
	}

	void pattern(const Matrix_IJV &other)
	{
		if constexpr(!A::always_equal) if(this->ator != other.ator) eslog::error("not implemented for unequal allocators");
		realloc(_allocated, other);
		_Matrix_IJV<T, I>::operator=(_allocated);
		this->rows = other.rows;
		this->cols = other.cols;
	}

	Matrix_Type type;
	Matrix_Shape shape;

protected:
	template <typename Type>
	void swap(Type &v, Type &u)
	{
		Type tmp = v; v = u; u = tmp;
	}

	void swap(_Matrix_IJV<T, I> &m, _Matrix_IJV<T, I> &n)
	{
		swap(m.nrows, n.nrows);
		swap(m.ncols, n.ncols);
		swap(m.nnz, n.nnz);
		swap(m.rows, n.rows);
		swap(m.cols, n.cols);
		swap(m.vals, n.vals);
	}

	void realloc(_Matrix_IJV<T, I> &m, I nrows, I ncols, I nnz)
	{
		if (m.nnz < nnz) {
			if (m.rows) { ator.deallocate(m.rows); m.rows = nullptr; }
			if (m.cols) { ator.deallocate(m.cols); m.cols = nullptr; }
			if (m.vals) { ator.deallocate(m.vals); m.vals = nullptr; }
			m.rows = ator.template allocate<I>(nnz);
			m.cols = ator.template allocate<I>(nnz);
			m.vals = ator.template allocate<T>(nnz);
		}
		m.nrows = nrows;
		m.ncols = ncols;
		m.nnz = nnz;
	}

	void realloc(_Matrix_IJV<T, I> &m, const _Matrix_IJV<T, I> &other)
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

	void clear(_Matrix_IJV<T, I> &m)
	{
		m.nrows = m.ncols = m.nnz = 0;
		if (m.rows) { ator.deallocate(m.rows); m.rows = nullptr; }
		if (m.cols) { ator.deallocate(m.cols); m.cols = nullptr; }
		if (m.vals) { ator.deallocate(m.vals); m.vals = nullptr; }
	}

	_Matrix_IJV<T, I> _allocated;
	A ator;
};

}

#endif /* SRC_MATH2_PRIMITIVES_MATRIX_IJV_H_ */
