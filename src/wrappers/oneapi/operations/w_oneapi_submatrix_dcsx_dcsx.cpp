
#ifdef HAVE_ONEAPI

#include "wrappers/oneapi/operations/w_oneapi_submatrix_dcsx_dcsx.h"

#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunused-local-typedef"
#pragma clang diagnostic ignored "-Wdeprecated-declarations"
#include <oneapi/dpl/execution>
#pragma clang diagnostic pop
#include <oneapi/dpl/async>
namespace onedpl = oneapi::dpl;

#include "wrappers/oneapi/common_oneapi_mgm.h"



namespace espreso {
namespace gpu {
namespace operations {



template<typename T, typename I>
struct functor_calc_start_end_ptrs_and_nnzperprim
{
    I * src_ptrs;
    I * src_idxs;
    I * start_ptrs;
    I * end_ptrs;
    I * nnz_per_prim;
    size_t start_primary;
    size_t start_secdary;
    size_t end_secdary;
    functor_calc_start_end_ptrs_and_nnzperprim(I * src_ptrs_, I * src_idxs_, I * start_ptrs_, I * end_ptrs_, I * nnz_per_prim_, size_t start_primary_, size_t start_secdary_, size_t end_secdary_)
        : src_ptrs(src_ptrs_)
        , src_idxs(src_idxs_)
        , start_ptrs(start_ptrs_)
        , end_ptrs(end_ptrs_)
        , nnz_per_prim(nnz_per_prim_)
        , start_primary(start_primary_)
        , start_secdary(start_secdary_)
        , end_secdary(end_secdary_)
    {}
    void operator()(sycl::nd_item<1> item) const {
        sycl::group g = item.get_group();
        I prim_dst = g.get_group_linear_range() - g.get_group_linear_id() - 1; // because my matrices are heavier towards the bottom, so assign the longest blocks with blockIdx.x=0, and hope that they get executed first
        I prim_src = prim_dst + start_primary;

        I start = src_ptrs[prim_src];
        I end = src_ptrs[prim_src+1];
        if(g.get_local_linear_id() == 0) {
            start_ptrs[prim_dst] = start;
            end_ptrs[prim_dst] = end;
        }
        group_barrier(g);
        for(I i = start + g.get_local_linear_id(); i < end-1; i += g.get_local_linear_range()) {
            I secdary_curr = src_idxs[i];
            I secdary_next = src_idxs[i+1];
            if(secdary_curr < start_secdary && start_secdary <= secdary_next) {
                start_ptrs[prim_dst] = i+1;
            }
            if(secdary_curr < end_secdary && end_secdary <= secdary_next) {
                end_ptrs[prim_dst] = i+1;
            }
        }
        if(g.get_local_linear_id() == 0) {
            group_barrier(g);
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
};



template<typename T, typename I>
struct functor_calc_dst_idxs_vals
{
    I * src_start_ptrs;
    I * src_end_ptrs;
    I * src_idxs;
    T * src_vals;
    I * dst_ptrs;
    I * dst_idxs;
    T * dst_vals;
    I start_primary;
    I start_secdary;
    functor_calc_dst_idxs_vals(I * src_start_ptrs_, I * src_end_ptrs_, I * src_idxs_, T * src_vals_, I * dst_ptrs_, I * dst_idxs_, T * dst_vals_, I start_primary_, I start_secdary_)
        : src_start_ptrs(src_start_ptrs_)
        , src_end_ptrs(src_end_ptrs_)
        , src_idxs(src_idxs_)
        , src_vals(src_vals_)
        , dst_ptrs(dst_ptrs_)
        , dst_idxs(dst_idxs_)
        , dst_vals(dst_vals_)
        , start_primary(start_primary_)
        , start_secdary(start_secdary_)
    {}
    void operator()(sycl::nd_item<1> item) const {
        sycl::group g = item.get_group();
        I prim_dst = g.get_group_linear_range() - g.get_group_linear_id() - 1; // because my matrices are heavier towards the bottom, so assign the longest blocks with blockIdx.x=0, and hope that they get executed first

        I src_start = src_start_ptrs[prim_dst];
        I src_end = src_end_ptrs[prim_dst];
        I dst_start = dst_ptrs[prim_dst];
        for(I i_src = src_start + g.get_local_linear_id(); i_src < src_end; i_src += g.get_local_linear_range()) {
            I offset = i_src - src_start;
            I i_dst = dst_start + offset;
            I iss = src_idxs[i_src];
            I isd = iss - start_secdary;
            T val = src_vals[i_src];
            dst_idxs[i_dst] = isd;
            dst_vals[i_dst] = val;
        }
    }
};



template<typename T, typename I>
void w_oneapi_submatrix_dcsx_dcsx<T,I>::internal_setup()
{
    size_t dst_size_primary = M_dst->get_size_primary();

    wss_pers_startptrs = utils::round_up(dst_size_primary * sizeof(I), gpu::mgm::get_natural_pitch_align());
    wss_pers_endptrs = utils::round_up(dst_size_primary * sizeof(I), gpu::mgm::get_natural_pitch_align());
    wss_pers_outptrs = utils::round_up((dst_size_primary + 1) * sizeof(I), gpu::mgm::get_natural_pitch_align());

    wss_internal = 0;
    wss_persistent = wss_pers_startptrs + wss_pers_endptrs + wss_pers_outptrs;
    wss_tmp_preprocess = wss_scan;
    wss_tmp_perform = 0;
}



template<typename T, typename I>
void w_oneapi_submatrix_dcsx_dcsx<T,I>::internal_preprocess(void * ws_tmp)
{
    size_t dst_size_primary = M_dst->get_size_primary();

    src_start_ptrs = reinterpret_cast<I*>(ws_persistent);
    src_end_ptrs = reinterpret_cast<I*>((char*)ws_persistent + wss_pers_startptrs);
    dst_ptrs = reinterpret_cast<I*>((char*)ws_persistent + wss_pers_startptrs + wss_pers_endptrs);

    int group_size = 256;
    sycl::nd_range<1> range(dst_size_primary * group_size, group_size);
    q->q.parallel_for(range, functor_calc_start_end_ptrs_and_nnzperprim<T,I>(M_src->ptrs, M_src->idxs, src_start_ptrs, src_end_ptrs, dst_ptrs, primary_start, secdary_start, secdary_end));

    onedpl::experimental::exclusive_scan_async(onedpl::execution::make_device_policy(q->q), dst_ptrs, dst_ptrs + dst_size_primary + 1, dst_ptrs, I{0}, std::plus<I>());
}



template<typename T, typename I>
void w_oneapi_submatrix_dcsx_dcsx<T,I>::internal_perform(void * /*ws_tmp*/)
{
    size_t dst_size_primary = M_dst->get_size_primary();

    q->q.template copy<I>(dst_ptrs, M_dst->ptrs, dst_size_primary + 1);

    int group_size = 256;
    sycl::nd_range<1> range(dst_size_primary * group_size, group_size);
    q->q.parallel_for(range, functor_calc_dst_idxs_vals<T,I>(src_start_ptrs, src_end_ptrs, M_src->idxs, M_src->vals, M_dst->ptrs, M_dst->idxs, M_dst->vals, primary_start, secdary_start));
}



#define INSTANTIATE_T_I(T,I) \
template class w_oneapi_submatrix_dcsx_dcsx<T,I>;

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
