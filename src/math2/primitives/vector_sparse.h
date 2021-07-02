
#ifndef SRC_MATH2_PRIMITIVES_VECTOR_SPARSE_H_
#define SRC_MATH2_PRIMITIVES_VECTOR_SPARSE_H_

namespace espreso {

struct _Vector_Sparse_Pattern {
	esint size, nnz, *indices;
};

template <typename T>
struct _Vector_Sparse_Vals {
	T *vals;
};

template <typename T>
struct _Vector_Sparse: public _Vector_Sparse_Pattern, public _Vector_Sparse_Vals<T> {

};

template <typename T>
class Vector_Sparse: public _Vector_Sparse<T>
{
public:
	Vector_Sparse(): _Vector_Sparse<T>{}, _allocated{}
	{

	}

	Vector_Sparse(const Vector_Sparse &other): _Vector_Sparse<T>{}, _allocated{}
	{
		realloc(_allocated, other.size, other.nnz);
		_Vector_Sparse<T>::operator=(_allocated);
		for (esint i = 0; i < other.nnz; ++i) {
			this->indices[i] = other.indices[i];
			this->vals[i] = other.vals[i];
		}
	}

	Vector_Sparse(Vector_Sparse &&other): _Vector_Sparse<T>{}, _allocated{}
	{
		swap(*this, other);
		swap(_allocated, other._allocated);
	}

	Vector_Sparse& operator=(const Vector_Sparse &other)
	{
		realloc(_allocated, other.size, other.nnz);
		_Vector_Sparse<T>::operator=(_allocated);
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
		_Vector_Sparse<T>::operator=(_allocated);
	}

	void resize(const Vector_Sparse &other)
	{
		resize(other.size, other.nnz);
	}

	void pattern(const Vector_Sparse &other)
	{
		realloc(_allocated, other);
		_Vector_Sparse<T>::operator=(_allocated);
		this->indices = other.indices;
	}

protected:
	template <typename Type>
	void swap(Type &v, Type &u)
	{
		Type tmp = v; v = u; u = tmp;
	}

	void swap(_Vector_Sparse<T> &v, _Vector_Sparse<T> &u)
	{
		swap(v.size, u.size);
		swap(v.nnz, u.nnz);
		swap(v.indices, u.indices);
		swap(v.vals, u.vals);
	}

	void realloc(_Vector_Sparse<T> &v, esint size, esint nnz)
	{
		if (v.nnz < nnz) {
			clear(v);
			v.indices = new esint[nnz];
			v.vals = new T[nnz];
		}
		v.size = size;
		v.nnz = nnz;
	}

	void realloc(_Vector_Sparse<T> &v, const _Vector_Sparse_Pattern &other)
	{
		if (v.indices) { delete[] v.indices; v.indices = nullptr; }

		if (v.nnz < other.nnz) {
			clear(v);
			v.vals = new T[other.nnz];
		}
		v.size = other.size;
		v.nnz = other.nnz;
	}

	void clear(_Vector_Sparse<T> &v)
	{
		if (v.indices) { delete[] v.indices; v.indices = nullptr; }
		if (v.vals) { delete[] v.vals; v.vals = nullptr; }
	}

	_Vector_Sparse<T> _allocated;
};

}

#endif /* SRC_MATH2_PRIMITIVES_VECTOR_SPARSE_H_ */
