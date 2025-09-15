
#ifdef HAVE_ONEAPI

#include "wrappers/oneapi/operations/w_oneapi_submatrix_ddnx_ddnx_noncontig.h"

#include "wrappers/oneapi/common_oneapi_mgm.h"
#include "wrappers/oneapi/common_oneblas.h"



namespace espreso {
namespace gpu {
namespace operations {



template<typename T, typename I, bool SUB_PRIMARY, bool SUB_SECDARY>
struct functor_do_submatrix
{
    I size_dst_primary;
    I size_dst_secdary;
    T * src;
    I ld_src;
    T * dst;
    I ld_dst;
    I * map_primary;
    I * map_secdary;
    functor_do_submatrix(I size_dst_primary_, I size_dst_secdary_, T * src_, I ld_src_, T * dst_, I ld_dst_, I * map_primary_, I * map_secdary_)
        : size_dst_primary(size_dst_primary_)
        , size_dst_secdary(size_dst_secdary_)
        , src(src_)
        , ld_src(ld_src_)
        , dst(dst_)
        , ld_dst(ld_dst_)
        , map_primary(map_primary_)
        , map_secdary(map_secdary_)
    {}
    void operator()(sycl::nd_item<2> item) const {
        sycl::group g = item.get_group();
        I ipd_start = g.get_group_id(0) * g.get_local_range(0) + g.get_local_id(0);
        I ipd_stride = g.get_local_range(0) * g.get_group_range(0);
        I isd_start = g.get_group_id(1) * g.get_local_range(1) + g.get_local_id(1);
        I isd_stride = g.get_local_range(1) * g.get_group_range(1);

        for(I ipd = ipd_start; ipd < size_dst_primary; ipd += ipd_stride) {
            I ips = ipd;
            if constexpr(SUB_PRIMARY) {
                ips = map_primary[ipd];
            }

            T * dst_prim = dst + ipd * ld_dst;
            T * src_prim = src + ips * ld_src;

            for(I isd = isd_start; isd < size_dst_secdary; isd += isd_stride) {
                I iss = isd;
                if constexpr(SUB_SECDARY) {
                    iss = map_secdary[isd];
                }
                dst_prim[isd] = src_prim[iss];
            }
        }
    }
};



template<typename T, typename I>
void w_oneapi_submatrix_ddnx_ddnx_noncontig<T,I>::internal_perform()
{
    VectorDenseView_new<I> * map_primary = ((d_M_src->order == 'R') ? d_row_map : d_col_map);
    VectorDenseView_new<I> * map_secdary = ((d_M_src->order == 'R') ? d_col_map : d_row_map);

    sycl::range local_range(8,32);
    sycl::range n(2,2);
    sycl::range global_range(utils::divide_round_up(d_M_dst->get_size_primary(), n[0]), utils::divide_round_up(d_M_dst->get_size_secdary(), n[1]));
    sycl::nd_range range(global_range, local_range);

    if(map_primary == nullptr && map_secdary == nullptr) {
        oneblas::row_major::omatcopy(q->q, onemkl::transpose::nontrans, d_M_src->get_size_primary(), d_M_src->get_size_secdary(), T{1}, d_M_src->vals, d_M_src->ld, d_M_dst->vals, d_M_dst->ld);
    }
    if(map_primary != nullptr && map_secdary == nullptr) {
        q->q.parallel_for(range, functor_do_submatrix<T,I,true,false>(d_M_dst->get_size_primary(), d_M_dst->get_size_secdary(), d_M_src->vals, d_M_src->ld, d_M_dst->vals, d_M_dst->ld, map_primary->vals, nullptr));
    }
    if(map_primary == nullptr && map_secdary != nullptr) {
        q->q.parallel_for(range, functor_do_submatrix<T,I,false,true>(d_M_dst->get_size_primary(), d_M_dst->get_size_secdary(), d_M_src->vals, d_M_src->ld, d_M_dst->vals, d_M_dst->ld, nullptr, map_secdary->vals));
    }
    if(map_primary != nullptr && map_secdary != nullptr) {
        q->q.parallel_for(range, functor_do_submatrix<T,I,true,true>(d_M_dst->get_size_primary(), d_M_dst->get_size_secdary(), d_M_src->vals, d_M_src->ld, d_M_dst->vals, d_M_dst->ld, map_primary->vals, map_secdary->vals));
    }
}



#define INSTANTIATE_T_I(T,I) \
template class w_oneapi_submatrix_ddnx_ddnx_noncontig<T,I>;

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
