
#ifdef HAVE_CUDA

#include "wrappers/cuda/operations/w_cuda_submatrix_dcsx_ddny.h"



namespace espreso {
namespace gpu {
namespace operations {



template<typename T, typename I>
__global__
static void submatrix_csx_dnx_vals(I * src_ptrs, I * src_idxs, T * src_vals, T * dst_vals, I dst_ld, I src_primary_start, I src_secdary_start, I src_secdary_end)
{
    I ips = blockIdx.x;
    I ipd = ips - src_primary_start;
    I start = src_ptrs[ips];
    I end = src_ptrs[ips+1];
    for(I i = start + threadIdx.x; i < end; i += blockDim.x) {
        I iss = src_idxs[i];
        if(iss >= src_secdary_end) {
            break;
        }
        T val = src_vals[i];
        I isd = iss - src_secdary_start;
        dst_vals[ipd * dst_ld + isd] = val;
    }
}



template<typename T, typename I>
__global__
static void submatrix_csx_dny_vals(I * src_ptrs, I * src_idxs, T * src_vals, T * dst_vals, I dst_ld, I src_primary_start, I src_secdary_start, I src_secdary_end)
{
    I ips = blockIdx.x;
    I isd = ips - src_primary_start;
    I start = src_ptrs[ips];
    I end = src_ptrs[ips+1];
    for(I i = start + threadIdx.x; i < end; i += blockDim.x) {
        I iss = src_idxs[i];
        if(iss >= src_secdary_end) {
            break;
        }
        T val = src_vals[i];
        I ipd = iss - src_secdary_start;
        dst_vals[ipd * dst_ld + isd] = val;
    }
}



template<typename T, typename I>
void w_cuda_submatrix_dcsx_ddny<T,I>::internal_setup()
{
    wss_internal = 0;
    wss_persistent = 0;
    wss_tmp_preprocess = 0;
    wss_tmp_perform = 0;
}



template<typename T, typename I>
void w_cuda_submatrix_dcsx_ddny<T,I>::internal_preprocess(void * /*ws_tmp*/)
{
    // no-op
}



template<typename T, typename I>
void w_cuda_submatrix_dcsx_ddny<T,I>::internal_perform(void * /*ws_tmp*/)
{
    CHECK(cudaMemset2D(M_dst->vals, M_dst->ld * sizeof(T), 0, M_dst->get_size_secdary() * sizeof(T), M_dst->get_size_primary()));

    if(M_src->order == M_dst->order) {
        submatrix_csx_dnx_vals<T,I><<<M_dst->get_size_primary(),256>>>(M_src->ptrs, M_src->idxs, M_src->vals, M_dst->vals, M_dst->ld, primary_start, secdary_start, secdary_end);
        CHECK(cudaPeekAtLastError());
    }
    else {
        submatrix_csx_dny_vals<T,I><<<M_dst->get_size_primary(),256>>>(M_src->ptrs, M_src->idxs, M_src->vals, M_dst->vals, M_dst->ld, primary_start, secdary_start, secdary_end);
        CHECK(cudaPeekAtLastError());
    }
}



#define INSTANTIATE_T_I(T,I) \
template class w_cuda_submatrix_dcsx_ddny<T,I>;

    #define INSTANTIATE_T(T) \
    INSTANTIATE_T_I(T,int32_t) \
    /* INSTANTIATE_T_I(T,int64_t) */

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
}
}

#endif
