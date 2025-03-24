
#ifdef HAVE_CUDA

#include "wrappers/cuda/operations/w_cuda_supermatrix_ddnx_ddnx_noncontig.h"

#include "wrappers/cuda/common_cuda_mgm.h"



namespace espreso {
namespace gpu {
namespace operations {



template<typename T, typename I>
__global__
static void do_supermatrix_primary(I size_src_primary, I size_secdary, T * src, I ld_src, T * dst, I ld_dst, I * map_primary)
{
    I ips = blockIdx.x;
    I ipd = map_primary[ips];

    T * dst_prim = dst + ipd * ld_dst;
    T * src_prim = src + ips * ld_src;

    for(I is = threadIdx.x; is < size_secdary; is += blockDim.x) {
        dst_prim[is] = src_prim[is];
    }
}



template<typename T, typename I>
__global__
static void do_supermatrix_secdary(I size_primary, I size_src_secdary, T * src, I ld_src, T * dst, I ld_dst, I * map_secdary)
{
    I ip = blockIdx.x;

    T * dst_prim = dst + ip * ld_dst;
    T * src_prim = src + ip * ld_src;

    for(I iss = threadIdx.x; iss < size_src_secdary; iss += blockDim.x) {
        I isd = map_secdary[iss];
        dst_prim[isd] = src_prim[iss];
    }
}



template<typename T, typename I>
__global__
static void do_supermatrix_both(I size_src_primary, I size_src_secdary, T * src, I ld_src, T * dst, I ld_dst, I * map_primary, I * map_secdary)
{
    I ips = blockIdx.x;
    I ipd = map_primary[ips];

    T * dst_prim = dst + ipd * ld_dst;
    T * src_prim = src + ips * ld_src;

    for(I iss = threadIdx.x; iss < size_src_secdary; iss += blockDim.x) {
        I isd = map_secdary[iss];
        dst_prim[isd] = src_prim[iss];
    }
}



template<typename T, typename I>
void w_cuda_supermatrix_ddnx_ddnx_noncontig<T,I>::internal_perform()
{
    VectorDenseView_new<I> * map_primary = ((d_M_src->order == 'R') ? d_row_map : d_col_map);
    VectorDenseView_new<I> * map_secdary = ((d_M_src->order == 'R') ? d_col_map : d_row_map);

    if(map_primary == nullptr && map_secdary == nullptr) {
        CHECK(cudaMemcpy2DAsync(d_M_dst->vals, d_M_dst->ld * sizeof(T), d_M_src->vals, d_M_src->ld * sizeof(T), d_M_src->get_size_secdary() * sizeof(T), d_M_src->get_size_primary(), cudaMemcpyDefault, q->stream));
    }
    if(map_primary != nullptr && map_secdary == nullptr) {
        int bpg = d_M_src->get_size_primary();
        int tpb = 256;
        do_supermatrix_primary<T,I><<<bpg,tpb,0,q->stream>>>(d_M_src->get_size_primary(), d_M_src->get_size_secdary(), d_M_src->vals, d_M_src->ld, d_M_dst->vals, d_M_dst->ld, map_primary->vals);
        CHECK(cudaPeekAtLastError());
    }
    if(map_primary == nullptr && map_secdary != nullptr) {
        int bpg = d_M_src->get_size_primary();
        int tpb = 256;
        do_supermatrix_secdary<T,I><<<bpg,tpb,0,q->stream>>>(d_M_src->get_size_primary(), d_M_src->get_size_secdary(), d_M_src->vals, d_M_src->ld, d_M_dst->vals, d_M_dst->ld, map_secdary->vals);
        CHECK(cudaPeekAtLastError());
    }
    if(map_primary != nullptr && map_secdary != nullptr) {
        int bpg = d_M_src->get_size_primary();
        int tpb = 256;
        do_supermatrix_both<T,I><<<bpg,tpb,0,q->stream>>>(d_M_src->get_size_primary(), d_M_src->get_size_secdary(), d_M_src->vals, d_M_src->ld, d_M_dst->vals, d_M_dst->ld, map_primary->vals, map_secdary->vals);
        CHECK(cudaPeekAtLastError());
    }
}



#define INSTANTIATE_T_I(T,I) \
template class w_cuda_supermatrix_ddnx_ddnx_noncontig<T,I>;

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
