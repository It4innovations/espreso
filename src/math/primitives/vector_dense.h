
#ifndef SRC_MATH2_PRIMITIVES_VECTOR_DENSE_H_
#define SRC_MATH2_PRIMITIVES_VECTOR_DENSE_H_

#include "basis/containers/allocators.h"
#include "esinfo/eslog.hpp"

namespace espreso {

template <typename T, typename I>
struct _Vector_Dense {
	I size;
	T *vals;
};

template <typename T, typename I = int, typename A = cpu_allocator>
class Vector_Dense: public _Vector_Dense<T, I>
{
public:
	Vector_Dense(const A &ator_ = A()): _Vector_Dense<T, I>{}, _allocated{}, ator(ator_)
	{

	}

	Vector_Dense(const Vector_Dense &other): _Vector_Dense<T, I>{}, _allocated{}, ator(other.ator)
	{
		realloc(_allocated, other.size);
		_Vector_Dense<T, I>::operator=(_allocated);
		for (I i = 0; i < other.size; ++i) {
			this->vals[i] = other.vals[i];
		}
	}

	Vector_Dense(Vector_Dense &&other): _Vector_Dense<T, I>{}, _allocated{}, ator(std::move(other.ator))
	{
		swap(*this, other);
		swap(_allocated, other._allocated);
	}

	Vector_Dense& operator=(const Vector_Dense &other) = delete;
//	{
//		realloc(_allocated, other.size);
//		_Vector_Dense<T, I>::operator=(_allocated);
//		for (esint i = 0; i < other.size; ++i) {
//			this->vals[i] = other.vals[i];
//		}
//		return *this;
//	}

	Vector_Dense& operator=(Vector_Dense &&other) = delete;
//	{
//		swap(*this, other);
//		swap(_allocated, other._allocated);
//		return *this;
//	}

	~Vector_Dense()
	{
		clear(_allocated);
	}

	void resize(I size)
	{
		realloc(_allocated, size);
		_Vector_Dense<T, I>::operator=(_allocated);
	}

	template<typename T2, typename I2, typename A2>
	void resize(const Vector_Dense<T2,I2,A2> &other)
	{
		resize(other.size);
	}

	template<typename T2>
	void pattern(const Vector_Dense<T2,I,A> &other)
	{
		if constexpr(!A::always_equal) if(this->ator != other.ator) eslog::error("not implemented for unequal allocators\n");
		realloc(_allocated, other.size);
		_Vector_Dense<T, I>::operator=(_allocated);
	}

	void shallowCopy(const Vector_Dense &other)
	{
		if constexpr(!A::always_equal) if(this->ator != other.ator) eslog::error("not implemented for unequal allocators\n");
		_Vector_Dense<T, I>::operator=(other);
	}

	void swap(Vector_Dense &other)
	{
		swap(*this, other);
		swap(_allocated, other._allocated);
		std::swap(this->ator, other.ator);
	}

	_Vector_Dense<T, I> _allocated;
	A ator;

protected:
	template <typename Type>
	void _swap(Type &v, Type &u)
	{
		Type tmp = v; v = u; u = tmp;
	}

	void swap(_Vector_Dense<T, I> &v, _Vector_Dense<T, I> &u)
	{
		_swap(v.size, u.size);
		_swap(v.vals, u.vals);
	}

	void realloc(_Vector_Dense<T, I> &v, I size)
	{
		if (v.size < size) {
			clear(v);
			v.vals = ator.template allocate<T>(size);
		}
		v.size = size;
	}

	void clear(_Vector_Dense<T, I> &v)
	{
		if (v.vals) { ator.deallocate(v.vals); v.vals = nullptr; }
	}
};

}

#endif /* SRC_MATH2_PRIMITIVES_VECTOR_DENSE_H_ */
