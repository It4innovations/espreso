
#ifndef SRC_MATH2_PRIMITIVES_MATRIX_CSR_H_
#define SRC_MATH2_PRIMITIVES_MATRIX_CSR_H_

#include "matrix_info.h"
#include "basis/containers/allocators.h"
#include "esinfo/eslog.hpp"

namespace espreso {

struct Matrix_CSR_Solver;

template <typename T, typename I>
struct _Matrix_CSR {
	I nrows, ncols, nnz, *rows, *cols;
	T *vals;
};

template <typename T, typename I = int, typename A = cpu_allocator>
class Matrix_CSR: public _Matrix_CSR<T, I>
{
public:
	Matrix_CSR(const A &ator_ = A()): _Matrix_CSR<T, I>{}, type{Matrix_Type::REAL_STRUCTURALLY_SYMMETRIC}, shape{Matrix_Shape::FULL}, _allocated{}, ator(ator_)
	{

	}

	Matrix_CSR(const Matrix_CSR &other): _Matrix_CSR<T, I>{}, type{Matrix_Type::REAL_STRUCTURALLY_SYMMETRIC}, shape{Matrix_Shape::FULL}, _allocated{}, ator(other.ator)
	{
		type = other.type;
		shape = other.shape;
		realloc(_allocated, other.nrows, other.ncols, other.nnz);
		_Matrix_CSR<T, I>::operator=(_allocated);
		for (I i = 0; i <= other.nrows; ++i) {
			this->rows[i] = other.rows[i];
		}
		for (I i = 0; i < other.nnz; ++i) {
			this->cols[i] = other.cols[i];
			this->vals[i] = other.vals[i];
		}
	}

	Matrix_CSR(Matrix_CSR &&other): _Matrix_CSR<T, I>{}, type{Matrix_Type::REAL_STRUCTURALLY_SYMMETRIC}, shape{Matrix_Shape::FULL}, _allocated{}, ator(std::move(other.ator))
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

	void resize(I nrows, I ncols, I nnz)
	{
		realloc(_allocated, nrows, ncols, nnz);
		_Matrix_CSR<T, I>::operator=(_allocated);
	}

	template<typename T2, typename I2, typename A2>
	void resize(const Matrix_CSR<T2,I2,A2> &other)
	{
		resize(other.nrows, other.ncols, other.nnz);
	}

	template<typename T2>
	void pattern(const Matrix_CSR<T2,I,A> &other)
	{
		if constexpr(!A::always_equal) if(this->ator != other.ator) eslog::error("not implemented for unequal allocators\n");
		realloc(_allocated, other);
		_Matrix_CSR<T, I>::operator=(_allocated);
		this->rows = other.rows;
		this->cols = other.cols;
	}

	void shallowCopy(const Matrix_CSR &other)
	{
		if constexpr(!A::always_equal) if(this->ator != other.ator) eslog::error("not implemented for unequal allocators\n");
		type = other.type;
		shape = other.shape;
		_Matrix_CSR<T, I>::operator=(other);
	}

    static size_t memoryRequirement(I nrows, I /*ncols*/, I nvals)
    {
        return (nrows+1) * sizeof(I) + nvals * sizeof(I) + nvals * sizeof(T);
    }

	Matrix_Type type;
	Matrix_Shape shape;
	_Matrix_CSR<T, I> _allocated;
	A ator;

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

	void realloc(_Matrix_CSR<T, I> &m, I nrows, I ncols, I nnz)
	{
		if (m.nrows < nrows) {
			if (m.rows) { ator.deallocate(m.rows); m.rows = nullptr; }
		}
		if (m.rows == nullptr) {
			m.rows = ator.template allocate<I>(nrows + 1);
		}
		if (m.nnz < nnz) {
			if (m.cols) { ator.deallocate(m.cols); m.cols = nullptr; }
			if (m.vals) { ator.deallocate(m.vals); m.vals = nullptr; }
		}
		if (m.cols == nullptr) {
			m.cols = ator.template allocate<I>(nnz);
		}
		if (m.vals == nullptr) {
			m.vals = ator.template allocate<T>(nnz);
		}
		m.nrows = nrows;
		m.ncols = ncols;
		m.nnz = nnz;
	}

	void realloc(_Matrix_CSR<T, I> &m, const _Matrix_CSR<T, I> &other)
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

	void clear(_Matrix_CSR<T, I> &m)
	{
		m.nrows = m.ncols = m.nnz = 0;
		if (m.rows) { ator.deallocate(m.rows); m.rows = nullptr; }
		if (m.cols) { ator.deallocate(m.cols); m.cols = nullptr; }
		if (m.vals) { ator.deallocate(m.vals); m.vals = nullptr; }
	}
};

}

#endif /* SRC_MATH2_PRIMITIVES_MATRIX_CSR_H_ */
