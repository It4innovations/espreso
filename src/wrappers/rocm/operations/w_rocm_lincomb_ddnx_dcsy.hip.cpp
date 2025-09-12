
#ifdef HAVE_ROCM

#include "wrappers/rocm/operations/w_rocm_lincomb_ddnx_dcsy.h"

#include "wrappers/rocm/common_rocm_mgm.h"
#include "wrappers/rocm/common_internal.hip.h"



namespace espreso {
namespace gpu {
namespace operations {



template<typename T, typename I>
__global__
static void w_rocm_lincomb_ddnx_dcsy_kernel_sameorder(I * sp_ptrs, I * sp_idxs, T * sp_vals, T * dn_vals, size_t dn_ld, utils::remove_complex_t<T> beta)
{
    size_t ip = blockIdx.x;

    T * dn_ptr = dn_vals + ip * dn_ld;

    I start = sp_ptrs[ip];
    I end = sp_ptrs[ip+1];
    for(I i = start + threadIdx.x; i < end; i += blockDim.x)
    {
        I is = sp_idxs[i];
        T val = sp_vals[i];
        myGenericAddScaled(&dn_ptr[is], val, beta);
    }
}



template<typename T, typename I>
__global__
static void w_rocm_lincomb_ddnx_dcsy_kernel_difforder(I * sp_ptrs, I * sp_idxs, T * sp_vals, T * dn_vals, size_t dn_ld, utils::remove_complex_t<T> beta)
{
    size_t ips = blockIdx.x;
    I isd = ips;

    I start = sp_ptrs[ips];
    I end = sp_ptrs[ips+1];
    for(I i = start + threadIdx.x; i < end; i += blockDim.x)
    {
        I iss = sp_idxs[i];
        T val = sp_vals[i];

        I ipd = iss;
        myGenericAddScaled(&dn_vals[ipd * dn_ld + isd], val, beta);
    }
}



template<typename T, typename I>
void w_rocm_lincomb_ddnx_dcsy<T,I>::internal_setup()
{
    wss_tmp_perform = 0;

    if(X != A || alpha != T{1}) eslog::error("only in-place supported now (A=X, alpha=1)\n");
}



template<typename T, typename I>
void w_rocm_lincomb_ddnx_dcsy<T,I>::internal_perform(void * /*ws_tmp*/)
{
    if(X != A || alpha != T{1}) eslog::error("only in-place supported now (A=X, alpha=1)\n");

    if(X->order == B->order) {
        w_rocm_lincomb_ddnx_dcsy_kernel_sameorder<T,I><<< B->get_size_primary(), 64, 0, q->stream >>>(B->ptrs, B->idxs, B->vals, X->vals, X->ld, beta);
        CHECK(hipPeekAtLastError());
    }
    else {
        w_rocm_lincomb_ddnx_dcsy_kernel_difforder<T,I><<< B->get_size_primary(), 64, 0, q->stream >>>(B->ptrs, B->idxs, B->vals, X->vals, X->ld, beta);
        CHECK(hipPeekAtLastError());
    }
}



#define INSTANTIATE_T_I(T,I) \
template class w_rocm_lincomb_ddnx_dcsy<T,I>;

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
