
#ifdef HAVE_ROCM

#include "wrappers/rocm/operations/w_rocm_submatrix_dcsx_ddny.h"

#include "wrappers/rocm/common_rocm_mgm.h"



namespace espreso {
namespace gpu {
namespace operations {

// TODO: do some preprocessing to find start and end indexes of each row
//       so I dont have to scan the whole row only to extract a section of it



template<typename T, typename I, bool SAME_ORDER>
__global__
static void submatrix_csx_dny_vals(I * src_ptrs, I * src_idxs, T * src_vals, T * dst_vals, I dst_ld, I src_primary_start, I src_secdary_start, I src_secdary_end)
{
    I ipd = blockIdx.x;
    I ips = ipd + src_primary_start;
    I start = src_ptrs[ips];
    I end = src_ptrs[ips+1];
    for(I i = start + threadIdx.x; i < end; i += blockDim.x) {
        I iss = src_idxs[i];
        if(iss < src_secdary_start) {
            continue;
        }
        if(iss >= src_secdary_end) {
            break;
        }
        I isd = iss - src_secdary_start;
        T val = src_vals[i];
        if constexpr(SAME_ORDER) {
            dst_vals[ipd * dst_ld + isd] = val;
        }
        else {
            dst_vals[isd * dst_ld + ipd] = val;
        }
    }
}



template<typename T, typename I>
void w_rocm_submatrix_dcsx_ddny<T,I>::internal_setup()
{
    wss_internal = 0;
    wss_persistent = 0;
    wss_tmp_preprocess = 0;
    wss_tmp_perform = 0;
}



template<typename T, typename I>
void w_rocm_submatrix_dcsx_ddny<T,I>::internal_preprocess(void * /*ws_tmp*/)
{
    // no-op
}



template<typename T, typename I>
void w_rocm_submatrix_dcsx_ddny<T,I>::internal_perform(void * /*ws_tmp*/)
{
    CHECK(hipMemset2D(M_dst->vals, M_dst->ld * sizeof(T), 0, M_dst->get_size_secdary() * sizeof(T), M_dst->get_size_primary()));

    if(M_src->order == M_dst->order) {
        submatrix_csx_dny_vals<T,I,true><<<M_dst->get_size_primary(),256,0,q->stream>>>(M_src->ptrs, M_src->idxs, M_src->vals, M_dst->vals, M_dst->ld, primary_start, secdary_start, secdary_end);
        CHECK(hipPeekAtLastError());
    }
    else {
        submatrix_csx_dny_vals<T,I,false><<<M_dst->get_size_primary(),256,0,q->stream>>>(M_src->ptrs, M_src->idxs, M_src->vals, M_dst->vals, M_dst->ld, primary_start, secdary_start, secdary_end);
        CHECK(hipPeekAtLastError());
    }
}



#define INSTANTIATE_T_I(T,I) \
template class w_rocm_submatrix_dcsx_ddny<T,I>;

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
