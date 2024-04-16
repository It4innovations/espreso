
#ifndef SRC_MATH2_PRIMITIVES_MATRIX_DENSE_H_
#define SRC_MATH2_PRIMITIVES_MATRIX_DENSE_H_

#include "matrix_info.h"
#include "slice.h"
#include "basis/containers/allocators.h"
#include "esinfo/eslog.hpp"

#include <limits>

namespace espreso {

template <typename T, typename I>
struct _Matrix_Dense {
    I nrows, ncols, nnz, maxnnz;
    T *vals;
};

template <typename T, typename I = int, typename A = cpu_allocator>
class Matrix_Dense: public _Matrix_Dense<T, I>
{
public:
    Matrix_Dense(const A &ator_ = A()): _Matrix_Dense<T, I>{}, type{Matrix_Type::REAL_STRUCTURALLY_SYMMETRIC}, shape{Matrix_Shape::FULL}, _allocated{}, ator(ator_)
    {

    }

    Matrix_Dense(const Matrix_Dense &other): _Matrix_Dense<T, I>{}, _allocated{}, ator(other.ator)
    {
        this->type = other.type;
        this->shape = other.shape;
        realloc(_allocated, this->shape, other.nrows, other.ncols);
        _Matrix_Dense<T, I>::operator=(_allocated);
        for (I i = 0; i < other.nnz; ++i) {
            this->vals[i] = other.vals[i];
        }
    }

    Matrix_Dense(Matrix_Dense &&other): _Matrix_Dense<T, I>{}, _allocated{}, ator(std::move(other.ator))
    {
        this->type = other.type;
        this->shape = other.shape;
        swap(*this, other);
        swap(_allocated, other._allocated);
    }

    Matrix_Dense& operator=(const Matrix_Dense &other) = delete;
//    {
//        this->type = other.type;
//        this->shape = other.shape;
//        realloc(_allocated, other.nrows, other.ncols);
//        _Matrix_Dense<T, I>::operator=(_allocated);
//        for (esint i = 0; i < other.nrows * other.ncols; ++i) {
//            this->vals[i] = other.vals[i];
//        }
//        return *this;
//    }

    Matrix_Dense& operator=(Matrix_Dense &&other) = delete;
//    {
//        this->type = other.type;
//        this->shape = other.shape;
//        swap(*this, other);
//        swap(_allocated, other._allocated);
//        return *this;
//    }

    ~Matrix_Dense()
    {
        clear(_allocated);
    }

    _Matrix_Dense<T, I>& allocated()
    {
        return _allocated;
    }

    void resize(I nrows, I ncols)
    {
        realloc(_allocated, shape, nrows, ncols);
        _Matrix_Dense<T, I>::operator=(_allocated);
    }

    void resize(I nrows_, I ncols_, I ld_)
    {
        if(ld_ < 0) ld_ = ncols_;
        resize(nrows_, ld_);
        this->ncols = ncols_;
    }

    template<typename T2, typename I2, typename A2>
    void resize(const Matrix_Dense<T2,I2,A2> &other)
    {
        resize(other.nrows, other.ncols, other.get_ld());
    }

    template<typename T2>
    void pattern(const Matrix_Dense<T2,I,A> &other)
    {
        if constexpr(!A::always_equal) if(this->ator != other.ator) eslog::error("not implemented for unequal allocators\n");
        realloc(_allocated, other.shape, other.nrows, other.ncols);
        _Matrix_Dense<T, I>::operator=(_allocated);
    }

    void slice(const Slice &rows, const Slice &cols)
    {
        submatrix[0] = rows;
        submatrix[1] = cols;
    }

    I get_ld() const
    {
        return _allocated.ncols;
    }

    static size_t memoryRequirement(I nrows, I ncols, I ld)
    {
        if(ld < 0) ld = ncols;
        return nrows * (ld * sizeof(T));
    }

    Matrix_Type type;
    Matrix_Shape shape;
    Slice submatrix[2];

protected:
    template <typename Type>
    void swap(Type &v, Type &u)
    {
        Type tmp = v; v = u; u = tmp;
    }

    void swap(_Matrix_Dense<T, I> &m, _Matrix_Dense<T, I> &n)
    {
        swap(m.nrows, n.nrows);
        swap(m.ncols, n.ncols);
        swap(m.vals, n.vals);
    }

    void realloc(_Matrix_Dense<T, I> &m, const Matrix_Shape &shape, I nrows, I ncols)
    {
        if((size_t)nrows * ncols > std::numeric_limits<I>::max()) eslog::error("matrix too large for the used integer type\n");
        I nnz;
        if (shape == Matrix_Shape::FULL) {
            nnz = nrows * ncols;
        } else {
            nnz = (nrows * ncols - nrows) / 2 + nrows;
        }
        if (m.maxnnz < nrows * ncols) {
            clear(m);
            m.vals = ator.template allocate<T>(nnz);
            m.maxnnz = nnz;
        }
        m.nrows = nrows;
        m.ncols = ncols;
        m.nnz = nnz;
    }

    void clear(_Matrix_Dense<T, I> &m)
    {
        m.nrows = m.ncols = 0;
        if (m.vals) { ator.deallocate(m.vals); m.vals = nullptr; }
        m.maxnnz = 0;
    }

    _Matrix_Dense<T, I> _allocated;
    A ator;
};


}

#endif /* SRC_MATH2_PRIMITIVES_MATRIX_DENSE_H_ */
