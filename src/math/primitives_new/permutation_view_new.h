
#ifndef SRC_MATH_PRIMITIVES_NEW_PERMUTATION_VIEW_NEW_H_
#define SRC_MATH_PRIMITIVES_NEW_PERMUTATION_VIEW_NEW_H_

#include "math/primitives_new/vector_base_new.h"
#include "math/primitives_new/allocator_new_base.h"
#include "math/primitives/permutation.h"



namespace espreso {



template<typename I>
struct PermutationView_new : public VectorBase_new
{
public: // the user promises not to modify these values (I don't want to implement getters everywhere)
    Allocator_new * ator = nullptr;
    I * dst_to_src = nullptr;   // idx_src = dst_to_src[idx_dst];   dst[i] = src[dst_to_src[i]];
    I * src_to_dst = nullptr;   // idx_dst = src_to_dst[idx_src];   dst[src_to_dst[i]] = src[i];
    bool was_set = false;
public:
    using VectorBase_new::size;
public:
    PermutationView_new() = default;
    PermutationView_new(const PermutationView_new &) = default;
    PermutationView_new(PermutationView_new &&) = default;
    PermutationView_new & operator=(const PermutationView_new &) = default;
    PermutationView_new & operator=(PermutationView_new &&) = default;
    virtual ~PermutationView_new() = default;
public:
    void set_view(size_t size_, I * dst_to_src_, I * src_to_dst_, Allocator_new * ator_)
    {
        if(was_set) eslog::error("can only set yet-uninitialized permutation view\n");
        size = size_;
        dst_to_src = dst_to_src_;
        src_to_dst = src_to_dst_;
        ator = ator_;
        was_set = true;
    }
public:
    PermutationView_new<I> get_inverse_view()
    {
        PermutationView_new<I> ret;
        ret.set_view(size, src_to_dst, dst_to_src, ator);
        return ret;
    }

    template<typename A>
    static PermutationView_new<I> from_old(const Permutation<I,A> & P_old)
    {
        PermutationView_new<I> P_new;
        P_new.set_view(P_old.size, P_old.dst_to_src, P_old.src_to_dst, AllocatorDummy_new::get_singleton(A::is_data_host_accessible, A::is_data_device_accessible));
        return P_new;
    }
    template<typename A>
    static Permutation<I,A> to_old(PermutationView_new<I> & P_new)
    {
        if(A::is_data_host_accessible != P_new.ator->is_data_accessible_cpu()) eslog::error("allocator access mismatch on cpu\n");
        if(A::is_data_device_accessible != P_new.ator->is_data_accessible_gpu()) eslog::error("allocator access mismatch on gpu\n");
        Permutation<I,A> P_old;
        P_old._allocated.size = P_new.size;
        P_old.size = P_new.size;
        P_old.dst_to_src = P_new.dst_to_src;
        P_old.src_to_dst = P_new.src_to_dst;
        return P_old;
    }

    static void invert(I * src, I * dst, size_t size)
    {
        for(size_t i = 0; i < size; i++) {
            dst[src[i]] = i;
        }
    }

    void print(const char * name = "")
    {
        if(!ator->is_data_accessible_cpu()) eslog::error("print is supported only for cpu-accessible matrices\n");
        eslog::info("Permutation %s, size %zu\n", name, size);
        eslog::info("dst_to_src: ");
        for(size_t i = 0; i < size; i++) {
            eslog::info(" %+11lld", (long long)dst_to_src[i]);
        }
        eslog::info("\n");
        eslog::info("src_to_dst: ");
        for(size_t i = 0; i < size; i++) {
            eslog::info(" %+11lld", (long long)src_to_dst[i]);
        }
        eslog::info("\n");
        fflush(stdout);
    }
};



}



#endif /* SRC_MATH_PRIMITIVES_NEW_PERMUTATION_VIEW_NEW_H_ */
