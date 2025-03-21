
#ifdef HAVE_CUDA

#include "wrappers/cuda/operations/w_cuda_permute_ddnx_ddnx.h"

#include "wrappers/cuda/common_cuda_mgm.h"



namespace espreso {
namespace gpu {
namespace operations {



template<typename T, typename I, bool PERM_PRIMARY, bool PERM_SECDARY>
__global__
static void permute_kernel(T * src, I ld_src, T * dst, I ld_dst, I size_secdary, I * perm_primary_dst_to_src, I * perm_secdary_dst_to_src)
{
    I ipd = blockIdx.x;
    I ips = ipd;
    if constexpr(PERM_PRIMARY) {
        ips = perm_primary_dst_to_src[ipd];
    }
    
    T * dstprim = dst + ipd * ld_dst;
    T * srcprim = src + ips * ld_src;

    for(I isd = threadIdx.x; isd < ipd; isd += blockDim.x) {
        I iss = isd;
        if constexpr(PERM_SECDARY) {
            iss = perm_secdary_dst_to_src[isd];
        }
        dstprim[isd] = srcprim[iss];
    }
}



template<typename T, typename I>
void w_cuda_permute_ddnx_ddnx<T,I>::internal_setup()
{
    wss_tmp_perform = 0;
}



template<typename T, typename I>
void w_cuda_permute_ddnx_ddnx<T,I>::internal_perform(void * ws_tmp)
{
    size_t size_primary = M_src->get_size_primary();
    size_t size_secdary = M_src->get_size_secdary();

    if(perm_primary == nullptr && perm_secdary == nullptr) {
        CHECK(cudaMemcpy2DAsync(M_dst->vals, M_dst->ld * sizeof(T), M_src->vals, M_src->ld * sizeof(T), M_src->get_size_secdary() * sizeof(T), M_src->get_size_primary(), cudaMemcpyDefault, q->stream));
    }
    if(perm_primary == nullptr && perm_secdary != nullptr) {
        permute_kernel<T,I,false,true ><<<size_primary,256,0,this->q->stream>>>(M_src->vals, M_src->ld, M_dst->vals, M_dst->ld, size_secdary, perm_primary->dst_to_src, perm_secdary->dst_to_src);
        CHECK(cudaPeekAtLastError());
    }
    if(perm_primary != nullptr && perm_secdary == nullptr) {
        permute_kernel<T,I,true, false><<<size_primary,256,0,this->q->stream>>>(M_src->vals, M_src->ld, M_dst->vals, M_dst->ld, size_secdary, perm_primary->dst_to_src, perm_secdary->dst_to_src);
        CHECK(cudaPeekAtLastError());
    }
    if(perm_primary != nullptr && perm_secdary != nullptr) {
        permute_kernel<T,I,true, true ><<<size_primary,256,0,this->q->stream>>>(M_src->vals, M_src->ld, M_dst->vals, M_dst->ld, size_secdary, perm_primary->dst_to_src, perm_secdary->dst_to_src);
        CHECK(cudaPeekAtLastError());
    }
}



#define INSTANTIATE_T_I(T,I) \
template class w_cuda_permute_ddnx_ddnx<T,I>;

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
