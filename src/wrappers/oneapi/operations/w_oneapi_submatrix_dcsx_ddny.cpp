
#ifdef HAVE_ONEAPI

#include "wrappers/oneapi/operations/w_oneapi_submatrix_dcsx_ddny.h"

#include "wrappers/oneapi/common_oneapi_mgm.h"
#include "wrappers/oneapi/common_oneblas.h"



namespace espreso {
namespace gpu {
namespace operations {

// TODO: do some preprocessing to find start and end indexes of each row
//       so I dont have to scan the whole row only to extract a section of it



template<typename T, typename I, bool SAME_ORDER>
struct functor_submatrix_csx_dny_vals
{
    I * src_ptrs;
    I * src_idxs;
    T * src_vals;
    T * dst_vals;
    I dst_ld;
    I src_primary_start;
    I src_secdary_start;
    I src_secdary_end;
    functor_submatrix_csx_dny_vals(I * src_ptrs_, I * src_idxs_, T * src_vals_, T * dst_vals_, I dst_ld_, I src_primary_start_, I src_secdary_start_, I src_secdary_end_)
        : src_ptrs(src_ptrs_)
        , src_idxs(src_idxs_)
        , src_vals(src_vals_)
        , dst_vals(dst_vals_)
        , dst_ld(dst_ld_)
        , src_primary_start(src_primary_start_)
        , src_secdary_start(src_secdary_start_)
        , src_secdary_end(src_secdary_end_)
    {}
    void operator()(sycl::nd_item<1> item) const {
        sycl::group g = item.get_group();
        I ipd = g.get_group_linear_id();
        I ips = ipd + src_primary_start;
        I start = src_ptrs[ips];
        I end = src_ptrs[ips+1];
        for(I i = start + g.get_local_linear_id(); i < end; i += g.get_local_linear_range()) {
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
};



template<typename T, typename I>
void w_oneapi_submatrix_dcsx_ddny<T,I>::internal_setup()
{
    wss_internal = 0;
    wss_persistent = 0;
    wss_tmp_preprocess = 0;
    wss_tmp_perform = 0;
}



template<typename T, typename I>
void w_oneapi_submatrix_dcsx_ddny<T,I>::internal_preprocess(void * /*ws_tmp*/)
{
    // no-op
}



template<typename T, typename I>
void w_oneapi_submatrix_dcsx_ddny<T,I>::internal_perform(void * /*ws_tmp*/)
{
    oneblas::row_major::imatcopy(q->q, onemkl::transpose::nontrans, M_dst->get_size_primary(), M_dst->get_size_secdary(), T{0}, M_dst->vals, M_dst->ld, M_dst->ld);

    int group_size = 256;
    sycl::nd_range<1> range(M_dst->get_size_primary() * group_size, group_size);

    if(M_src->order == M_dst->order) {
        q->q.parallel_for(range, functor_submatrix_csx_dny_vals<T,I,true>(M_src->ptrs, M_src->idxs, M_src->vals, M_dst->vals, M_dst->ld, primary_start, secdary_start, secdary_end));
    }
    else {
        q->q.parallel_for(range, functor_submatrix_csx_dny_vals<T,I,false>(M_src->ptrs, M_src->idxs, M_src->vals, M_dst->vals, M_dst->ld, primary_start, secdary_start, secdary_end));
    }
}



#define INSTANTIATE_T_I(T,I) \
template class w_oneapi_submatrix_dcsx_ddny<T,I>;

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
