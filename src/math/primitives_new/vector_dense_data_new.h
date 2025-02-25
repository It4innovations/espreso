
#ifndef SRC_MATH_PRIMITIVES_NEW_VECTOR_DENSE_NEW_H_
#define SRC_MATH_PRIMITIVES_NEW_VECTOR_DENSE_NEW_H_

#include "math/primitives_new/vector_base_view_new.h"
#include "math/primitives_new/allocator_new.h"



namespace espreso {



template<typename T>
struct VectorDenseData_new : public VectorDenseView_new<T>
{
public: // the user promises not to modify these values (I don't want to implement getters everywhere)
    Allocator_new * ator = nullptr;
public:
    using VectorDenseView_new::vals;
    using VectorBase_new::size;
public:
    VectorDenseData_new() = default;
    VectorDenseData_new(const VectorDenseData_new &) = delete;
    VectorDenseData_new(VectorDenseData_new && other)
    {
        std::swap(static_cast<VectorDenseView_new&>(*this), static_cast<VectorDenseView_new>(other));
        std::swap(ator, other.ator);
    }
    VectorDenseData_new & operator=(const VectorDenseData_new &) = delete;
    VectorDenseData_new & operator=(VectorDenseData_new && other)
    {
        if(this == &other) return;
        std::swap(static_cast<VectorDenseView_new&>(*this), static_cast<VectorDenseView_new>(other));
        std::swap(ator, other.ator);
        other.free();
    }
    virtual ~VectorDenseData_new()
    {
        free();
    }
public:
    void set(size_t size_, Allocator_new * ator_)
    {
        if(was_set) eslog::error("can only set yet-uninitialized vector data\n");
        size = size_;
        ator = ator_;
        was_set = true;
    }
    void alloc()
    {
        if(!was_set) eslog::error("vector has not been set\n");
        if(ator == nullptr) eslog::error("vector has not been set\n");
        if(vals != nullptr) eslog::error("vector already contains data\n");
        if(size > 0) {
            vals = ator->alloc<T>(size);
        }
    }
    void free()
    {
        if(ator != nullptr) {
            ator->free(vals);
        }
    }
    void clear()
    {
        free();
        size = 0;
        ator = nullptr;
        was_set = false;
    }
};



}



#endif /* SRC_MATH_PRIMITIVES_NEW_VECTOR_DENSE_NEW_H_ */
