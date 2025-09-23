
#ifdef HAVE_CUDA

#include "wrappers/cuda/operations/w_cuda_submatrix_ddnx_ddnx_noncontig.h"

#include "wrappers/cuda/common_cuda_mgm.h"



namespace espreso {
namespace gpu {
namespace operations {



template<typename T, typename I>
__global__
static void do_submatrix_primary(size_t size_dst_primary, size_t size_secdary, T * src, size_t ld_src, T * dst, size_t ld_dst, I * map_primary)
{
    I ipd_start = blockIdx.y * blockDim.y + threadIdx.y;
    I ipd_stride = blockDim.y * gridDim.y;
    I is_start = blockIdx.x * blockDim.x + threadIdx.x;
    I is_stride = blockDim.x * gridDim.x;

    for(I ipd = ipd_start; ipd < size_dst_primary; ipd += ipd_stride) {
        I ips = map_primary[ipd];

        T * dst_prim = dst + ipd * ld_dst;
        T * src_prim = src + ips * ld_src;

        for(I is = is_start; is < size_secdary; is += is_stride) {
            dst_prim[is] = src_prim[is];
        }
    }
}



template<typename T, typename I>
__global__
static void do_submatrix_secdary(size_t size_primary, size_t size_dst_secdary, T * src, size_t ld_src, T * dst, size_t ld_dst, I * map_secdary)
{
    I ip_start = blockIdx.y * blockDim.y + threadIdx.y;
    I ip_stride = blockDim.y * gridDim.y;
    I isd_start = blockIdx.x * blockDim.x + threadIdx.x;
    I isd_stride = blockDim.x * gridDim.x;

    for(I ip = ip_start; ip < size_primary; ip += ip_stride) {
        T * dst_prim = dst + ip * ld_dst;
        T * src_prim = src + ip * ld_src;

        for(I isd = isd_start; isd < size_dst_secdary; isd += isd_stride) {
            I iss = map_secdary[isd];
            dst_prim[isd] = src_prim[iss];
        }
    }
}



template<typename T, typename I>
__global__
static void do_submatrix_both(size_t size_dst_primary, size_t size_dst_secdary, T * src, size_t ld_src, T * dst, size_t ld_dst, I * map_primary, I * map_secdary)
{
    I ipd_start = blockIdx.y * blockDim.y + threadIdx.y;
    I ipd_stride = blockDim.y * gridDim.y;
    I isd_start = blockIdx.x * blockDim.x + threadIdx.x;
    I isd_stride = blockDim.x * gridDim.x;

    for(I ipd = ipd_start; ipd < size_dst_primary; ipd += ipd_stride) {
        I ips = map_primary[ipd];

        T * dst_prim = dst + ipd * ld_dst;
        T * src_prim = src + ips * ld_src;

        for(I isd = isd_start; isd < size_dst_secdary; isd += isd_stride) {
            I iss = map_secdary[isd];
            dst_prim[isd] = src_prim[iss];
        }
    }
}



template<typename T, typename I>
void w_cuda_submatrix_ddnx_ddnx_noncontig<T,I>::internal_perform()
{
    if(d_M_dst->nrows == 0 || d_M_dst->ncols == 0) {
        return;
    }

    VectorDenseView_new<I> * map_primary = ((d_M_src->order == 'R') ? d_row_map : d_col_map);
    VectorDenseView_new<I> * map_secdary = ((d_M_src->order == 'R') ? d_col_map : d_row_map);

    dim3 tpb(32,8);
    dim3 n(2,2);
    dim3 bpg((d_M_dst->get_size_secdary() - 1) / (tpb.x * n.x) + 1, (d_M_dst->get_size_primary() - 1) / (tpb.y * n.y) + 1);

    if(map_primary == nullptr && map_secdary == nullptr) {
        CHECK(cudaMemcpy2DAsync(d_M_dst->vals, d_M_dst->ld * sizeof(T), d_M_src->vals, d_M_src->ld * sizeof(T), d_M_src->get_size_secdary() * sizeof(T), d_M_src->get_size_primary(), cudaMemcpyDefault, q->stream));
    }
    if(map_primary != nullptr && map_secdary == nullptr) {
        do_submatrix_primary<T,I><<<bpg,tpb,0,q->stream>>>(d_M_dst->get_size_primary(), d_M_dst->get_size_secdary(), d_M_src->vals, d_M_src->ld, d_M_dst->vals, d_M_dst->ld, map_primary->vals);
        CHECK(cudaPeekAtLastError());
    }
    if(map_primary == nullptr && map_secdary != nullptr) {
        do_submatrix_secdary<T,I><<<bpg,tpb,0,q->stream>>>(d_M_dst->get_size_primary(), d_M_dst->get_size_secdary(), d_M_src->vals, d_M_src->ld, d_M_dst->vals, d_M_dst->ld, map_secdary->vals);
        CHECK(cudaPeekAtLastError());
    }
    if(map_primary != nullptr && map_secdary != nullptr) {
        do_submatrix_both<T,I><<<bpg,tpb,0,q->stream>>>(d_M_dst->get_size_primary(), d_M_dst->get_size_secdary(), d_M_src->vals, d_M_src->ld, d_M_dst->vals, d_M_dst->ld, map_primary->vals, map_secdary->vals);
        CHECK(cudaPeekAtLastError());
    }
}



#define INSTANTIATE_T_I(T,I) \
template class w_cuda_submatrix_ddnx_ddnx_noncontig<T,I>;

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
