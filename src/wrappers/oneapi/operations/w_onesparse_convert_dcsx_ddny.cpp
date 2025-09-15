
#ifdef HAVE_ONEAPI

#include "wrappers/oneapi/operations/w_onesparse_convert_dcsx_ddny.h"

#include "wrappers/oneapi/common_oneblas.h"



namespace espreso {
namespace gpu {
namespace operations {



template<typename T, typename I, char order>
struct kernel_functor_sp2dn
{
    I * sparse_rowptrs;
    I * sparse_colidxs;
    T * sparse_vals;
    T * dense_vals;
    I dense_ld;
    kernel_functor_sp2dn(I*sr, I*sc, T*sv, T*dv, I ld)
        : sparse_rowptrs(sr)
        , sparse_colidxs(sc)
        , sparse_vals(sv)
        , dense_vals(dv)
        , dense_ld(ld)
    {}
    void operator()(sycl::nd_item<1> item) const {
        sycl::group g = item.get_group();
        I row = g.get_group_linear_id();
        I start = sparse_rowptrs[row];
        I end = sparse_rowptrs[row+1];
        for(I i = start + g.get_local_linear_id(); i < end; i += g.get_local_linear_range()) {
            I col = sparse_colidxs[i];
            T val = sparse_vals[i];
            if constexpr(order == 'R') dense_vals[row * dense_ld + col] = val;
            if constexpr(order == 'C') dense_vals[row + dense_ld * col] = val;
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
    {
        const I * src_ptrs = M_src->ptrs;
        const I * src_idxs = M_src->idxs;
        const T * src_vals = M_src->vals;
        T * dst_vals = M_dst->vals;
        size_t dst_ld = M_dst->ld;
        q->q.parallel_for(range, [=](sycl::nd_item<1> item) {
            sycl::group g = item.get_group();
            I ip = g.get_group_linear_id();
            T * dst_ptr = dst_vals + ip * dst_ld;
            I start = src_ptrs[ip];
            I end = src_ptrs[ip+1];
            for(I i = start + g.get_local_linear_id(); i < end; i += g.get_local_linear_range()) {
                I is = src_idxs[i];
                T val = src_vals[i];
                dst_ptr[is] = val;
            }
        });
    }
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
