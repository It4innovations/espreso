
#ifdef HAVE_ONEAPI

#include "wrappers/oneapi/operations/w_onesparse_convert_dcsx_ddny.h"

#include "wrappers/oneapi/common_oneblas.h"



namespace espreso {
namespace gpu {
namespace operations {



template<typename T, typename I, bool SAME_ORDER>
struct kernel_functor_sp2dn
{
    I * sparse_ptrs;
    I * sparse_idxs;
    T * sparse_vals;
    T * dense_vals;
    I dense_ld;
    kernel_functor_sp2dn(I*sp, I*si, T*sv, T*dv, I ld)
        : sparse_ptrs(sp)
        , sparse_idxs(si)
        , sparse_vals(sv)
        , dense_vals(dv)
        , dense_ld(ld)
    {}
    void operator()(sycl::nd_item<1> item) const {
        sycl::group g = item.get_group();
        I ip = g.get_group_linear_id();
        I start = sparse_ptrs[ip];
        I end = sparse_ptrs[ip+1];
        for(I i = start + g.get_local_linear_id(); i < end; i += g.get_local_linear_range()) {
            I is = sparse_idxs[i];
            T val = sparse_vals[i];
            if constexpr( SAME_ORDER) dense_vals[ip * dense_ld + is] = val;
            if constexpr(!SAME_ORDER) dense_vals[is * dense_ld + ip] = val;
        }
    }
};



template<typename T, typename I>
w_onesparse_convert_dcsx_ddny<T,I>::w_onesparse_convert_dcsx_ddny() = default;



template<typename T, typename I>
w_onesparse_convert_dcsx_ddny<T,I>::~w_onesparse_convert_dcsx_ddny() = default;



template<typename T, typename I>
void w_onesparse_convert_dcsx_ddny<T,I>::internal_setup()
{
    wss_internal = 0;
    wss_persistent = 0;
    wss_tmp_preprocess = 0;
    wss_tmp_perform = 0;
}



template<typename T, typename I>
void w_onesparse_convert_dcsx_ddny<T,I>::internal_preprocess(void * /*ws_tmp*/)
{
    // no-op
}



template<typename T, typename I>
void w_onesparse_convert_dcsx_ddny<T,I>::internal_perform(void * ws_tmp)
{
    oneblas::row_major::imatcopy(q->q, onemkl::transpose::nontrans, M_dst->get_size_primary(), M_dst->get_size_secdary(), T{0}, M_dst->vals, M_dst->ld, M_dst->ld);

    int workitems_in_workgroup = 64;
    sycl::nd_range<1> range(sycl::range<1>(M_src->get_size_primary() * workitems_in_workgroup), sycl::range<1>(workitems_in_workgroup));

    if(M_src->order == M_dst->order) q->q.parallel_for(range, kernel_functor_sp2dn<T,I,true>(M_src->ptrs, M_src->idxs, M_src->vals, M_dst->vals, M_dst->ld));
    if(M_src->order != M_dst->order) q->q.parallel_for(range, kernel_functor_sp2dn<T,I,false>(M_src->ptrs, M_src->idxs, M_src->vals, M_dst->vals, M_dst->ld));
}



#define INSTANTIATE_T_I(T,I) \
template class w_onesparse_convert_dcsx_ddny<T,I>;

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
