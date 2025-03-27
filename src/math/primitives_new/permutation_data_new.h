
#ifndef SRC_MATH_PRIMITIVES_NEW_PERMUTATION_DATA_NEW_H_
#define SRC_MATH_PRIMITIVES_NEW_PERMUTATION_DATA_NEW_H_

#include "math/primitives_new/permutation_view_new.h"
#include "math/primitives_new/allocator_new.h"



namespace espreso {



template<typename T>
struct PermutationData_new : public PermutationView_new<T>
{
public: // the user promises not to modify these values (I don't want to implement getters everywhere)
    Allocator_new * ator = nullptr;
public:
    using PermutationView_new<T>::dst_to_src;
    using PermutationView_new<T>::src_to_dst;
    using PermutationView_new<T>::was_set;
    using VectorBase_new::size;
public:
    PermutationData_new() = default;
    PermutationData_new(const PermutationData_new &) = delete;
    PermutationData_new(PermutationData_new && other)
    {
        std::swap(*static_cast<PermutationView_new<T>*>(this), *static_cast<PermutationView_new<T>*>(&other));
        std::swap(ator, other.ator);
    }
    PermutationData_new & operator=(const PermutationData_new &) = delete;
    PermutationData_new & operator=(PermutationData_new && other)
    {
        if(this != &other) {
            std::swap(*static_cast<PermutationView_new<T>*>(this), *static_cast<PermutationView_new<T>*>(&other));
            std::swap(ator, other.ator);
            other.free();
        }
        return *this;
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
public:
    size_t get_memory_impact()
    {
        size_t mem = size * sizeof(T);
        mem = utils::round_up(mem, ator->get_align());
        size_t total_mem = 2 * mem;
        return total_mem;
    }
};



}



#endif /* SRC_MATH_PRIMITIVES_NEW_PERMUTATION_DATA_NEW_H_ */
