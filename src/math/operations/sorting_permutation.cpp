
#include "math/operations/sorting_permutation.h"

#include "math/primitives_new/allocator_new.h"
#include "math/primitives_new/vector_dense_data_new.h"
#include "basis/utilities/stacktimer.h"



namespace espreso {
namespace math {
namespace operations {



template<typename T, typename I>
void sorting_permutation<T,I>::set_vector(VectorDenseView_new<T> * vec_)
{
    if(vec != nullptr) eslog::error("vector is already set\n");

    vec = vec_;
}



template<typename T, typename I>
void sorting_permutation<T,I>::set_permutation(PermutationView_new<I> * perm_)
{
    if(perm != nullptr) eslog::error("permutation is already set\n");

    perm = perm_;
}



template<typename T, typename I>
void sorting_permutation<T,I>::perform()
{
    stacktimer::push("sorting_permutation::perform");

    if(vec == nullptr) eslog::error("vector is not set\n");
    if(perm == nullptr) eslog::error("permutation is not set\n");
    if(!vec->ator->is_data_accessible_cpu()) eslog::error("vector must be cpu-accessible\n");
    if(!perm->ator->is_data_accessible_cpu()) eslog::error("permutation must be cpu-accessible\n");
    if(vec->size != perm->size) eslog::error("vector and permutation sizes dont match\n");

    struct idx_val { I idx; T val; };

    VectorDenseData_new<idx_val> idxsvals;
    idxsvals.set(vec->size, AllocatorCPU_new::get_singleton());
    idxsvals.alloc();

    I size = vec->size;
    for(I i = 0; i < size; i++) {
        idxsvals.vals[i].idx = i;
        idxsvals.vals[i].val = vec->vals[i];
    }

    std::sort(idxsvals.vals, idxsvals.vals + idxsvals.size, [](const idx_val & l, const idx_val & r){ return l.val < r.val; });

    for(I i = 0; i < size; i++) {
        perm->dst_to_src[i] = idxsvals.vals[i].idx;
    }
    PermutationView_new<I>::invert(perm->dst_to_src, perm->src_to_dst, perm->size);

    idxsvals.clear();

    stacktimer::pop();
}



template<typename T, typename I>
void sorting_permutation<T,I>::do_all(VectorDenseView_new<T> * vec, PermutationView_new<I> * perm)
{
    sorting_permutation<T,I> instance;
    instance.set_vector(vec);
    instance.set_permutation(perm);
    instance.perform();
}



#define INSTANTIATE_T_I(T,I) \
template class sorting_permutation<T,I>;

    #define INSTANTIATE_T(T) \
    INSTANTIATE_T_I(T,int32_t) \
    /* INSTANTIATE_T_I(T,int64_t) */ \
    INSTANTIATE_T_I(T,size_t)

        #define INSTANTIATE \
        INSTANTIATE_T(int32_t) \
        /* INSTANTIATE_T(int64_t) */ \
        /* INSTANTIATE_T(size_t) */ \
        /* INSTANTIATE_T(float) */ \
        /* INSTANTIATE_T(double) */

            INSTANTIATE

        #undef INSTANTIATE
    #undef INSTANTIATE_T
#undef INSTANTIATE_T_I



}
}
}

