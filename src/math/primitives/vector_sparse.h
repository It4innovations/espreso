
#ifndef SRC_MATH2_PRIMITIVES_VECTOR_SPARSE_H_
#define SRC_MATH2_PRIMITIVES_VECTOR_SPARSE_H_

#include "basis/containers/allocators.h"

namespace espreso {

template <typename T, typename I>
struct _Vector_Sparse {
	I size, nnz, *indices;
	T *vals;
};

template <typename T, typename I = int, template<typename> typename A = cpu_allocator>
class Vector_Sparse: public _Vector_Sparse<T, I>
{
public:
	Vector_Sparse(): _Vector_Sparse<T, I>{}, touched(false), _allocated{}
	{

	}

	Vector_Sparse(const Vector_Sparse &other): _Vector_Sparse<T, I>{}, touched(false), _allocated{}
	{
		realloc(_allocated, other.size, other.nnz);
		_Vector_Sparse<T, I>::operator=(_allocated);
		for (esint i = 0; i < other.nnz; ++i) {
			this->indices[i] = other.indices[i];
			this->vals[i] = other.vals[i];
		}
	}

	Vector_Sparse(Vector_Sparse &&other): _Vector_Sparse<T, I>{}, touched(false), _allocated{}
	{
		swap(*this, other);
		swap(_allocated, other._allocated);
	}

	Vector_Sparse& operator=(const Vector_Sparse &other)
	{
		realloc(_allocated, other.size, other.nnz);
		_Vector_Sparse<T, I>::operator=(_allocated);
		for (esint i = 0; i < other.nnz; ++i) {
			this->indices[i] = other.indices[i];
			this->vals[i] = other.vals[i];
		}
		return *this;
	}

	Vector_Sparse& operator=(Vector_Sparse &&other)
	{
		swap(*this, other);
		swap(_allocated, other._allocated);
		return *this;
	}

	~Vector_Sparse()
	{
		clear(_allocated);
	}

	void resize(esint size, esint nnz)
	{
		realloc(_allocated, size, nnz);
		_Vector_Sparse<T, I>::operator=(_allocated);
	}

	void resize(const Vector_Sparse &other)
	{
		resize(other.size, other.nnz);
	}

	void pattern(const Vector_Sparse &other)
	{
		realloc(_allocated, other);
		_Vector_Sparse<T, I>::operator=(_allocated);
		this->indices = other.indices;
	}

	bool touched;
	_Vector_Sparse<T, I> _allocated;

protected:
	template <typename Type>
	void swap(Type &v, Type &u)
	{
		Type tmp = v; v = u; u = tmp;
	}

	void swap(_Vector_Sparse<T, I> &v, _Vector_Sparse<T, I> &u)
	{
		swap(v.size, u.size);
		swap(v.nnz, u.nnz);
		swap(v.indices, u.indices);
		swap(v.vals, u.vals);
	}

	void realloc(_Vector_Sparse<T, I> &v, esint size, esint nnz)
	{
		if (v.nnz < nnz) {
			clear(v);
			v.indices = new esint[nnz];
			v.vals = new T[nnz];
		}
		v.size = size;
		v.nnz = nnz;
	}

	void realloc(_Vector_Sparse<T, I> &v, const _Vector_Sparse<T, I> &other)
	{
		if (v.indices) { delete[] v.indices; v.indices = nullptr; }

		if (v.nnz < other.nnz) {
			clear(v);
			v.vals = new T[other.nnz];
		}
		v.size = other.size;
		v.nnz = other.nnz;
	}

	void clear(_Vector_Sparse<T, I> &v)
	{
		if (v.indices) { delete[] v.indices; v.indices = nullptr; }
		if (v.vals) { delete[] v.vals; v.vals = nullptr; }
	}
};

}

#endif /* SRC_MATH2_PRIMITIVES_VECTOR_SPARSE_H_ */
