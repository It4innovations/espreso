
#ifdef HAVE_ROCM

#include "wrappers/rocm/operations/w_rocm_submatrix_dcsx_dcsx.h"

#include "wrappers/rocm/common_rocm_mgm.h"

#include <rocprim/rocprim.hpp>



namespace espreso {
namespace gpu {
namespace operations {



template<typename T, typename I>
__global__
static void calc_start_end_ptrs_and_nnzperprim(I * src_ptrs, I * src_idxs, I * start_ptrs, I * end_ptrs, I * nnz_per_prim, size_t start_primary, size_t start_secdary, size_t end_secdary)
{
    I prim_dst = gridDim.x - blockIdx.x - 1; // because my matrices are heavier towards the bottom, so assign the longest blocks with blockIdx.x=0, and hope that they get executed first
    I prim_src = prim_dst + start_primary;

    I start = src_ptrs[prim_src];
    I end = src_ptrs[prim_src+1];
    if(threadIdx.x == 0) {
        start_ptrs[prim_dst] = start;
        end_ptrs[prim_dst] = end;
    }
    __syncthreads();
    for(I i = start + threadIdx.x; i < end-1; i += blockDim.x) {
        I secdary_curr = src_idxs[i];
        I secdary_next = src_idxs[i+1];
        if(secdary_curr < start_secdary && start_secdary <= secdary_next) {
            start_ptrs[prim_dst] = i+1;
        }
        if(secdary_curr < end_secdary && end_secdary <= secdary_next) {
            end_ptrs[prim_dst] = i+1;
        }
    }
    __syncthreads();
    if(threadIdx.x == 0) {
        if(start != end) {
            if(end_secdary <= src_idxs[start]) {
                end_ptrs[prim_dst] = start;
            }
            if(start_secdary > src_idxs[end-1]) {
                start_ptrs[prim_dst] = end;
            }
        }
        nnz_per_prim[prim_dst] = end_ptrs[prim_dst] - start_ptrs[prim_dst];
    }
}



template<typename T, typename I>
__global__
static void calc_dst_idxs_vals(I * src_start_ptrs, I * src_end_ptrs, I * src_idxs, T * src_vals, I * dst_ptrs, I * dst_idxs, T * dst_vals, I start_primary, I start_secdary)
{
    I prim_dst = gridDim.x - blockIdx.x - 1; // because my matrices are heavier towards the bottom, so assign the longest blocks with blockIdx.x=0, and hope that they get executed first

    I src_start = src_start_ptrs[prim_dst];
    I src_end = src_end_ptrs[prim_dst];
    I dst_start = dst_ptrs[prim_dst];
    for(I i_src = src_start + threadIdx.x; i_src < src_end; i_src += blockDim.x) {
        I offset = i_src - src_start;
        I i_dst = dst_start + offset;
        I iss = src_idxs[i_src];
        I isd = iss - start_secdary;
        T val = src_vals[i_src];
        dst_idxs[i_dst] = isd;
        dst_vals[i_dst] = val;
    }
}



template<typename T, typename I>
void w_rocm_submatrix_dcsx_dcsx<T,I>::internal_setup()
{
    size_t dst_size_primary = M_dst->get_size_primary();

    wss_pers_startptrs = utils::round_up(dst_size_primary * sizeof(I), gpu::mgm::get_natural_pitch_align());
    wss_pers_endptrs = utils::round_up(dst_size_primary * sizeof(I), gpu::mgm::get_natural_pitch_align());
    wss_pers_outptrs = utils::round_up((dst_size_primary + 1) * sizeof(I), gpu::mgm::get_natural_pitch_align());

    CHECK(rocprim::exclusive_scan(nullptr, wss_scan, (I*)nullptr, (I*)nullptr, I{0}, dst_size_primary + 1, [] __device__ (I a, I b){return a+b;}, q->stream));

    wss_internal = 0;
    wss_persistent = wss_pers_startptrs + wss_pers_endptrs + wss_pers_outptrs;
    wss_tmp_preprocess = wss_scan;
    wss_tmp_perform = 0;
}



template<typename T, typename I>
void w_rocm_submatrix_dcsx_dcsx<T,I>::internal_preprocess(void * ws_tmp)
{
    size_t dst_size_primary = M_dst->get_size_primary();

    src_start_ptrs = reinterpret_cast<I*>(ws_persistent);
    src_end_ptrs = reinterpret_cast<I*>((char*)ws_persistent + wss_pers_startptrs);
    dst_ptrs = reinterpret_cast<I*>((char*)ws_persistent + wss_pers_startptrs + wss_pers_endptrs);

    calc_start_end_ptrs_and_nnzperprim<T,I><<<dst_size_primary,256,0,q->stream>>>(M_src->ptrs, M_src->idxs, src_start_ptrs, src_end_ptrs, dst_ptrs, primary_start, secdary_start, secdary_end);
    CHECK(hipPeekAtLastError());

    CHECK(rocprim::exclusive_scan(ws_tmp, wss_scan, dst_ptrs, dst_ptrs, I{0}, dst_size_primary + 1, [] __device__ (I a, I b){return a+b;}, q->stream));
}



template<typename T, typename I>
void w_rocm_submatrix_dcsx_dcsx<T,I>::internal_perform(void * /*ws_tmp*/)
{
    size_t dst_size_primary = M_dst->get_size_primary();

    CHECK(hipMemcpyAsync(M_dst->ptrs, dst_ptrs, (dst_size_primary + 1) * sizeof(I), hipMemcpyDefault, q->stream));

    calc_dst_idxs_vals<T,I><<<dst_size_primary,256,0,q->stream>>>(src_start_ptrs, src_end_ptrs, M_src->idxs, M_src->vals, M_dst->ptrs, M_dst->idxs, M_dst->vals, primary_start, secdary_start);
    CHECK(hipPeekAtLastError());
}



#define INSTANTIATE_T_I(T,I) \
template class w_rocm_submatrix_dcsx_dcsx<T,I>;

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
