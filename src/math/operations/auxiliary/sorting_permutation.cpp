
#include "math/operations/auxiliary/sorting_permutation.h"

#include "math/primitives_new/vector_dense_data_new.h"



namespace espreso {
namespace math {
namespace operations {



template<typename T, typename I>
void sorting_permutation<T,I>::set_vector(VectorDenseView_new<T> * vec_)
{
    vec = vec_;
}



template<typename T, typename I>
void sorting_permutation<T,I>::set_permutation(PermutationView_new<I> * perm_)
{
    perm = perm_;
}



template<typename T, typename I>
void sorting_permutation<T,I>::perform()
{
    if(vec == nullptr) eslog::error("vector is not set\n");
    if(perm == nullptr) eslog::error("permutation is not set\n");
    if(vec->size != perm->size) eslog::error("vector and permutation sizes dont match\n");

    struct idx_val { I idx; T val; };

    VectorDenseData_new<idx_val> idxsvals;
    idxsvals.set(vec->size, AllocatorCPU_new::get_singleton());
    idxsvals.alloc();

    I size = vec->size;
    for(I i = 0; i < size; i++) {
        idxsvals[i].idx = i;
        idxsvals[i].val = vec->vals[i];
    }

    std::sort(idxsvals->vals, idxsvals->vals + idxsvals->size, [](const idx_val & l, const idx_val & r){ return l.val < r.val; });

    for(I i = 0; i < size; i++) {
        perm->dst_to_src[i] = idxsvals->vals[i].idx;
    }
    PermutationView_new<I>::invert(perm->dst_to_src, perm->src_to_dst, perm->size);

    idxsvals.clear();
}



template<typename T, typename I>
void sorting_permutation<T,I>::do_all(VectorDenseView_new<T> * vec, PermutationView_new<I> * perm)
{
    sorting_permutation<T,I> instance;
    instance.set_vector(vec);
    instance.set_permutation(perm);
    instance.perform();
}



}
}
}

