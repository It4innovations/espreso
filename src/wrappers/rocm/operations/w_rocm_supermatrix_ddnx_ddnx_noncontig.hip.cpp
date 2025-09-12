
#ifdef HAVE_ROCM

#include "wrappers/rocm/operations/w_rocm_supermatrix_ddnx_ddnx_noncontig.h"

#include "wrappers/rocm/common_rocm_mgm.h"
#include "wrappers/rocm/common_internal.hip.h"



namespace espreso {
namespace gpu {
namespace operations {



template<typename T, typename I, typename OP>
__global__
static void do_supermatrix_none(size_t size_primary, size_t size_secdary, T * __restrict__ src, size_t ld_src, T * __restrict__ dst, size_t ld_dst, OP op)
{
    I ip_start = blockIdx.y * blockDim.y + threadIdx.y;
    I ip_stride = blockDim.y * gridDim.y;
    I is_start = blockIdx.x * blockDim.x + threadIdx.x;
    I is_stride = blockDim.x * gridDim.x;

    for(I ip = ip_start; ip < size_primary; ip += ip_stride) {
        T * dst_prim = dst + ip * ld_dst;
        T * src_prim = src + ip * ld_src;

        for(I is = is_start; is < size_secdary; is += is_stride) {
            op(dst_prim[is], src_prim[is]);
        }
    }
}



template<typename T, typename I, typename OP>
__global__
static void do_supermatrix_primary(size_t size_src_primary, size_t size_secdary, T * __restrict__ src, size_t ld_src, T * __restrict__ dst, size_t ld_dst, I * __restrict__ map_primary, OP op)
{
    I ips_start = blockIdx.y * blockDim.y + threadIdx.y;
    I ips_stride = blockDim.y * gridDim.y;
    I is_start = blockIdx.x * blockDim.x + threadIdx.x;
    I is_stride = blockDim.x * gridDim.x;

    for(I ips = ips_start; ips < size_src_primary; ips += ips_stride) {
        I ipd = map_primary[ips];

        T * dst_prim = dst + ipd * ld_dst;
        T * src_prim = src + ips * ld_src;

        for(I is = is_start; is < size_secdary; is += is_stride) {
            op(dst_prim[is], src_prim[is]);
        }
    }
}



template<typename T, typename I, typename OP>
__global__
static void do_supermatrix_secdary(size_t size_primary, size_t size_src_secdary, T * __restrict__ src, size_t ld_src, T * __restrict__ dst, size_t ld_dst, I * __restrict__ map_secdary, OP op)
{
    I ip_start = blockIdx.y * blockDim.y + threadIdx.y;
    I ip_stride = blockDim.y * gridDim.y;
    I iss_start = blockIdx.x * blockDim.x + threadIdx.x;
    I iss_stride = blockDim.x * gridDim.x;

    for(I ip = ip_start; ip < size_primary; ip += ip_stride) {
        T * dst_prim = dst + ip * ld_dst;
        T * src_prim = src + ip * ld_src;

        for(I iss = iss_start; iss < size_src_secdary; iss += iss_stride) {
            I isd = map_secdary[iss];
            op(dst_prim[isd], src_prim[iss]);
        }
    }
}



template<typename T, typename I, typename OP>
__global__
static void do_supermatrix_both(size_t size_src_primary, size_t size_src_secdary, T * __restrict__ src, size_t ld_src, T * __restrict__ dst, size_t ld_dst, I * __restrict__ map_primary, I * __restrict__ map_secdary, OP op)
{
    I ips_start = blockIdx.y * blockDim.y + threadIdx.y;
    I ips_stride = blockDim.y * gridDim.y;
    I iss_start = blockIdx.x * blockDim.x + threadIdx.x;
    I iss_stride = blockDim.x * gridDim.x;

    for(I ips = ips_start; ips < size_src_primary; ips += ips_stride) {
        I ipd = map_primary[ips];

        T * dst_prim = dst + ipd * ld_dst;
        T * src_prim = src + ips * ld_src;

        for(I iss = iss_start; iss < size_src_secdary; iss += iss_stride) {
            I isd = map_secdary[iss];
            op(dst_prim[isd], src_prim[iss]);
        }
    }
}



template<typename T, typename I, typename OP>
static void do_supermatrix_launch_kernel_2(hipStream_t stream, size_t size_src_primary, size_t size_src_secdary, T * src, size_t ld_src, T * dst, size_t ld_dst, I * map_primary, I * map_secdary, OP op)
{
    dim3 tpb(64,4);
    dim3 n(2,2);
    dim3 bpg((size_src_secdary - 1) / (tpb.x * n.x) + 1, (size_src_primary - 1) / (tpb.y * n.y) + 1);

    if(map_primary == nullptr && map_secdary == nullptr) {
        do_supermatrix_none<T,I><<<bpg,tpb,0,stream>>>(size_src_primary, size_src_secdary, src, ld_src, dst, ld_dst, op);
        CHECK(hipPeekAtLastError());
    }
    if(map_primary != nullptr && map_secdary == nullptr) {
        do_supermatrix_primary<T,I><<<bpg,tpb,0,stream>>>(size_src_primary, size_src_secdary, src, ld_src, dst, ld_dst, map_primary, op);
        CHECK(hipPeekAtLastError());
    }
    if(map_primary == nullptr && map_secdary != nullptr) {
        do_supermatrix_secdary<T,I><<<bpg,tpb,0,stream>>>(size_src_primary, size_src_secdary, src, ld_src, dst, ld_dst, map_secdary, op);
        CHECK(hipPeekAtLastError());
    }
    if(map_primary != nullptr && map_secdary != nullptr) {
        do_supermatrix_both<T,I><<<bpg,tpb,0,stream>>>(size_src_primary, size_src_secdary, src, ld_src, dst, ld_dst, map_primary, map_secdary, op);
        CHECK(hipPeekAtLastError());
    }
}



template<typename T, typename I>
static void do_supermatrix_launch_kernel_1(hipStream_t stream, size_t size_src_primary, size_t size_src_secdary, T * src, size_t ld_src, T * dst, size_t ld_dst, I * map_primary, I * map_secdary, typename supermatrix_ddnx_ddnx_noncontig<T,I>::mode mode_val)
{
    using mode = typename supermatrix_ddnx_ddnx_noncontig<T,I>::mode;

    switch(mode_val) {
        case mode::assign: {
            auto op = [] __device__ (T & dst, T src){ dst = src; };
            do_supermatrix_launch_kernel_2(stream, size_src_primary, size_src_secdary, src, ld_src, dst, ld_dst, map_primary, map_secdary, op);
            break;
        }
        case mode::accumulate: {
            auto op = [] __device__ (T & dst, T src){ myGenericAdd(&dst, src); };
            do_supermatrix_launch_kernel_2(stream, size_src_primary, size_src_secdary, src, ld_src, dst, ld_dst, map_primary, map_secdary, op);
            break;
        }
        case mode::accumulate_atomic: {
            auto op = [] __device__ (T & dst, T src){ myAtomicAdd(&dst, src); };
            do_supermatrix_launch_kernel_2(stream, size_src_primary, size_src_secdary, src, ld_src, dst, ld_dst, map_primary, map_secdary, op);
            break;
        }
        default: {
            eslog::error("wrong mode\n");
        }
    }
}



template<typename T, typename I>
void w_rocm_supermatrix_ddnx_ddnx_noncontig<T,I>::internal_perform()
{
    VectorDenseView_new<I> * map_primary = ((d_M_src->order == 'R') ? d_row_map : d_col_map);
    VectorDenseView_new<I> * map_secdary = ((d_M_src->order == 'R') ? d_col_map : d_row_map);
    I * map_primary_vals = ((map_primary != nullptr) ? map_primary->vals : nullptr);
    I * map_secdary_vals = ((map_secdary != nullptr) ? map_secdary->vals : nullptr);

    do_supermatrix_launch_kernel_1<T,I>(q->stream, d_M_src->get_size_primary(), d_M_src->get_size_secdary(), d_M_src->vals, d_M_src->ld, d_M_dst->vals, d_M_dst->ld, map_primary_vals, map_secdary_vals, mode_val);
}



#define INSTANTIATE_T_I(T,I) \
template class w_rocm_supermatrix_ddnx_ddnx_noncontig<T,I>;

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
