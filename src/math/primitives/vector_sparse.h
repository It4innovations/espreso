
#ifndef SRC_MATH2_PRIMITIVES_VECTOR_SPARSE_H_
#define SRC_MATH2_PRIMITIVES_VECTOR_SPARSE_H_

#include "basis/containers/allocators.h"
#include "esinfo/eslog.hpp"

namespace espreso {

template <typename T, typename I>
struct _Vector_Sparse {
    I size, nnz, *indices;
    T *vals;
};

template <typename T, typename I = int, typename A = cpu_allocator>
class Vector_Sparse: public _Vector_Sparse<T, I>
{
public:
    Vector_Sparse(const A &ator_ = A()): _Vector_Sparse<T, I>{}, touched(false), _allocated{}, ator(ator_)
    {

    }

    Vector_Sparse(const Vector_Sparse &other): _Vector_Sparse<T, I>{}, touched(false), _allocated{}, ator(other.ator)
    {
        realloc(_allocated, other.size, other.nnz);
        _Vector_Sparse<T, I>::operator=(_allocated);
        for (I i = 0; i < other.nnz; ++i) {
            this->indices[i] = other.indices[i];
            this->vals[i] = other.vals[i];
        }
    }

    Vector_Sparse(Vector_Sparse &&other): _Vector_Sparse<T, I>{}, touched(false), _allocated{}, ator(std::move(other.ator))
    {
        swap(*this, other);
        swap(_allocated, other._allocated);
    }

    Vector_Sparse& operator=(const Vector_Sparse &other)
    {
        if constexpr(!A::always_equal) if(this->ator != other.ator) eslog::error("not implemented for unequal allocators\n");
        realloc(_allocated, other.size, other.nnz);
        _Vector_Sparse<T, I>::operator=(_allocated);
        for (I i = 0; i < other.nnz; ++i) {
            this->indices[i] = other.indices[i];
            this->vals[i] = other.vals[i];
        }
        return *this;
    }

    Vector_Sparse& operator=(Vector_Sparse &&other)
    {
        swap(*this, other);
        swap(_allocated, other._allocated);
        std::swap(this->ator, other.ator);
        return *this;
    }

    ~Vector_Sparse()
    {
        clear(_allocated);
    }

    void resize(I size, I nnz)
    {
        realloc(_allocated, size, nnz);
        _Vector_Sparse<T, I>::operator=(_allocated);
    }

    template<typename T2, typename I2, typename A2>
    void resize(const Vector_Sparse<T2,I2,A2> &other)
    {
        resize(other.size, other.nnz);
    }

    template<typename T2>
    void pattern(const Vector_Sparse<T2,I,A> &other)
    {
        if constexpr(!A::always_equal) if(this->ator != other.ator) eslog::error("not implemented for unequal allocators\n");
        realloc(_allocated, other);
        _Vector_Sparse<T, I>::operator=(_allocated);
        this->indices = other.indices;
    }

    void shallowCopy(const Vector_Sparse &other)
    {
        if constexpr(!A::always_equal) if(this->ator != other.ator) eslog::error("not implemented for unequal allocators\n");
        _Vector_Sparse<T, I>::operator=(other);
    }

    bool touched;
    _Vector_Sparse<T, I> _allocated;
    A ator;

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

    void realloc(_Vector_Sparse<T, I> &v, I size, I nnz)
    {
        if (v.nnz < nnz) {
            clear(v);
            v.indices = ator.template allocate<I>(nnz);
            v.vals = ator.template allocate<T>(nnz);
        }
        v.size = size;
        v.nnz = nnz;
    }

    void realloc(_Vector_Sparse<T, I> &v, const _Vector_Sparse<T, I> &other)
    {
        if (v.indices) { ator.deallocate(v.indices); v.indices = nullptr; }

        if (v.nnz < other.nnz) {
            clear(v);
            v.vals = ator.template allocate<T>(other.nnz);
        }
        v.size = other.size;
        v.nnz = other.nnz;
    }

    void clear(_Vector_Sparse<T, I> &v)
    {
        if (v.indices) { ator.deallocate(v.indices); v.indices = nullptr; }
        if (v.vals) { ator.deallocate(v.vals); v.vals = nullptr; }
    }
};

}

#endif /* SRC_MATH2_PRIMITIVES_VECTOR_SPARSE_H_ */
