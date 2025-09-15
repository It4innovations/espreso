
#ifdef HAVE_ONEAPI

#include "wrappers/oneapi/operations/w_oneapi_permute_ddnx_ddnx.h"

#include "wrappers/oneapi/common_oneapi_mgm.h"
#include "wrappers/oneapi/common_oneblas.h"



namespace espreso {
namespace gpu {
namespace operations {



template<typename T, typename I, bool PERM_PRIMARY, bool PERM_SECDARY>
struct permute_functor
{
    T * src;
    T * dst;
    I * perm_primary_dst_to_src;
    I * perm_secdary_dst_to_src;
    size_t ld_dst;
    size_t ld_src;
    size_t size_secdary;
    permute_functor(T * src_, size_t ld_src_, T * dst_, size_t ld_dst_, size_t size_secdary_, I * perm_primary_dst_to_src_, I * perm_secdary_dst_to_src_)
        : src(src_)
        , dst(dst_)
        , perm_primary_dst_to_src(perm_primary_dst_to_src_)
        , perm_secdary_dst_to_src(perm_secdary_dst_to_src_)
        , ld_dst(ld_dst_)
        , ld_src(ld_src_)
        , size_secdary(size_secdary_)
    {}
    void operator()(sycl::nd_item<1> item) const {
        sycl::group g = item.get_group();
        I ipd = g.get_group_linear_id();
        I ips = ipd;
        if constexpr(PERM_PRIMARY) {
            ips = perm_primary_dst_to_src[ipd];
        }
        
        T * dstprim = dst + ipd * ld_dst;
        T * srcprim = src + ips * ld_src;

        for(I isd = g.get_local_linear_id(); isd < size_secdary; isd += g.get_local_linear_range()) {
            I iss = isd;
            if constexpr(PERM_SECDARY) {
                iss = perm_secdary_dst_to_src[isd];
            }
            dstprim[isd] = srcprim[iss];
        }
    }
};



template<typename T, typename I>
void w_oneapi_permute_ddnx_ddnx<T,I>::internal_setup()
{
    wss_tmp_perform = 0;
}



template<typename T, typename I>
void w_oneapi_permute_ddnx_ddnx<T,I>::internal_perform(void * ws_tmp)
{
    size_t size_primary = M_src->get_size_primary();
    size_t size_secdary = M_src->get_size_secdary();
    int workitems_in_workgroup = 256;
    sycl::nd_range<1> range(sycl::range<1>(size_primary * workitems_in_workgroup), sycl::range<1>(workitems_in_workgroup));

    if(perm_primary == nullptr && perm_secdary == nullptr) {
        oneblas::row_major::omatcopy(q->q, onemkl::transpose::nontrans, size_primary, size_secdary, T{1}, M_src->vals, M_src->ld, M_dst->vals, M_dst->ld);
    }
    if(perm_primary == nullptr && perm_secdary != nullptr) {
        q->q.parallel_for(range, permute_functor<T,I,false,true>(M_src->vals, M_src->ld, M_dst->vals, M_dst->ld, size_secdary, perm_primary->dst_to_src, perm_secdary->dst_to_src));
    }
    if(perm_primary != nullptr && perm_secdary == nullptr) {
        q->q.parallel_for(range, permute_functor<T,I,true,false>(M_src->vals, M_src->ld, M_dst->vals, M_dst->ld, size_secdary, perm_primary->dst_to_src, perm_secdary->dst_to_src));
    }
    if(perm_primary != nullptr && perm_secdary != nullptr) {
        q->q.parallel_for(range, permute_functor<T,I,true,true>(M_src->vals, M_src->ld, M_dst->vals, M_dst->ld, size_secdary, perm_primary->dst_to_src, perm_secdary->dst_to_src));
    }
}



#define INSTANTIATE_T_I(T,I) \
template class w_oneapi_permute_ddnx_ddnx<T,I>;

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
