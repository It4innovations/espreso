
#ifndef SRC_MATH_PRIMITIVES_NEW_VECTOR_DENSE_NEW_H_
#define SRC_MATH_PRIMITIVES_NEW_VECTOR_DENSE_NEW_H_

#include "math/primitives_new/permutation_view_new.h"
#include "math/primitives_new/allocator_new.h"



namespace espreso {



template<typename T>
struct PermutationData_new : public PermutationView_new<T>
{
public: // the user promises not to modify these values (I don't want to implement getters everywhere)
    Allocator_new * ator = nullptr;
public:
    using PermutationView_new::dst_to_src;
    using PermutationView_new::src_to_dst;
    using VectorBase_new::size;
public:
    PermutationData_new() = default;
    PermutationData_new(const PermutationData_new &) = delete;
    PermutationData_new(PermutationData_new && other)
    {
        std::swap(static_cast<PermutationView_new&>(*this), static_cast<PermutationView_new>(other));
        std::swap(ator, other.ator);
    }
    PermutationData_new & operator=(const PermutationData_new &) = delete;
    PermutationData_new & operator=(PermutationData_new && other)
    {
        if(this == &other) return;
        std::swap(static_cast<PermutationView_new&>(*this), static_cast<PermutationView_new>(other));
        std::swap(ator, other.ator);
        other.free();
    }
    virtual ~PermutationData_new()
    {
        free();
    }
public:
    void set(size_t size_, Allocator_new * ator_)
    {
        if(was_set) eslog::error("can only set yet-uninitialized permutation data\n");
        size = size_;
        ator = ator_;
        was_set = true;
    }
    void alloc()
    {
        if(!was_set) eslog::error("permutation has not been set\n");
        if(ator == nullptr) eslog::error("permutation has not been set\n");
        if(dst_to_src != nullptr || src_to_dst != nullptr) eslog::error("permutation already contains data\n");
        if(size > 0) {
            dst_to_src = ator->alloc<T>(size);
            src_to_dst = ator->alloc<T>(size);
        }
    }
    void free()
    {
        if(ator != nullptr) {
            ator->free(dst_to_src);
            ator->free(src_to_dst);
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
