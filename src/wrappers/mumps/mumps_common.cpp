
#include "wrappers/mumps/mumps_common.h"

#include "math/primitives_new/allocator_new.h"



namespace espreso {



template<typename T, typename I>
void mumps_helper_csx_to_ijv(MatrixCsxView_new<T,I> & M, VectorDenseData_new<I> & ijv_rowidxs, VectorDenseData_new<I> & ijv_colidxs)
{
    VectorDenseData_new<I> idxs_primary;
    VectorDenseData_new<I> idxs_secdary;
    idxs_primary.set(M.nnz, AllocatorCPU_new::get_singleton());
    idxs_secdary.set(M.nnz, AllocatorCPU_new::get_singleton());
    idxs_primary.alloc();
    idxs_secdary.alloc();

    for(size_t ip = 0; ip < M.get_size_primary(); ip++) {
        I start = M.ptrs[ip];
        I end = M.ptrs[ip+1];
        for(I i = start; i < end; i++) {
            idxs_primary.vals[i] = ip + 1;
        }
    }

    for(size_t i = 0; i < M.nnz; i++) {
        idxs_secdary.vals[i] = M.idxs[i] + 1;
    }

    if(M.order == 'R') {
        ijv_rowidxs = std::move(idxs_primary);
        ijv_colidxs = std::move(idxs_secdary);
    }
    if(M.order == 'C') {
        ijv_rowidxs = std::move(idxs_secdary);
        ijv_colidxs = std::move(idxs_primary);
    }
}



#define INSTANTIATE_T_I(T,I) \
template void mumps_helper_csx_to_ijv<T,I>(MatrixCsxView_new<T,I> & M, VectorDenseData_new<I> & ijv_rowidxs, VectorDenseData_new<I> & ijv_colidxs);

    #define INSTANTIATE_T(T) \
    INSTANTIATE_T_I(T, int32_t) \
    /* INSTANTIATE_T_I(T, int64_t) */

        #define INSTANTIATE \
        /* INSTANTIATE_T(float) */ \
        INSTANTIATE_T(double) \
        /* INSTANTIATE_T(std::complex<float>) */ \
        INSTANTIATE_T(std::complex<double>)

            INSTANTIATE

        #undef INSTANTIATE
    #undef INSTANTIATE_T
#undef INSTANTIATE_T_I



}

