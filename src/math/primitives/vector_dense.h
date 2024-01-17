
#ifndef SRC_MATH2_PRIMITIVES_VECTOR_DENSE_H_
#define SRC_MATH2_PRIMITIVES_VECTOR_DENSE_H_

#include "basis/containers/allocators.h"

namespace espreso {

template <typename T, typename I>
struct _Vector_Dense {
	I size;
	T *vals;
};

template <typename T, typename I = int, template<typename> typename A = cpu_allocator>
class Vector_Dense: public _Vector_Dense<T, I>
{
public:
	Vector_Dense(): _Vector_Dense<T, I>{}, _allocated{}
	{

	}

	Vector_Dense(const Vector_Dense &other): _Vector_Dense<T, I>{}, _allocated{}
	{
		realloc(_allocated, other.size);
		_Vector_Dense<T, I>::operator=(_allocated);
		for (esint i = 0; i < other.size; ++i) {
			this->vals[i] = other.vals[i];
		}
	}

	Vector_Dense(Vector_Dense &&other): _Vector_Dense<T, I>{}, _allocated{}
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

	void resize(esint size)
	{
		realloc(_allocated, size);
		_Vector_Dense<T, I>::operator=(_allocated);
	}

	void resize(const Vector_Dense &other)
	{
		resize(other.size);
	}

	void pattern(const Vector_Dense &other)
	{
		realloc(_allocated, other.size);
		_Vector_Dense<T, I>::operator=(_allocated);
	}

	void swap(Vector_Dense &other)
	{
		swap(*this, other);
		swap(_allocated, other._allocated);
	}

	_Vector_Dense<T, I> _allocated;

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

	void realloc(_Vector_Dense<T, I> &v, esint size)
	{
		if (v.size < size) {
			clear(v);
			v.vals = new T[size];
		}
		v.size = size;
	}

	void clear(_Vector_Dense<T, I> &v)
	{
		if (v.vals) { delete[] v.vals; v.vals = nullptr; }
	}
};

}

#endif /* SRC_MATH2_PRIMITIVES_VECTOR_DENSE_H_ */
