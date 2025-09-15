
#ifdef HAVE_ONEAPI

#include "wrappers/oneapi/operations/w_oneapi_lincomb_ddnx_dcsy.h"

#include "wrappers/oneapi/common_oneapi_mgm.h"



namespace espreso {
namespace gpu {
namespace operations {



template<typename T, typename I>
void w_oneapi_lincomb_ddnx_dcsy<T,I>::internal_setup()
{
    wss_tmp_perform = 0;

    if(X != A || alpha != T{1}) eslog::error("only in-place supported now (A=X, alpha=1)\n");
}



template<typename T, typename I>
void w_oneapi_lincomb_ddnx_dcsy<T,I>::internal_perform(void * /*ws_tmp*/)
{
    if(X != A || alpha != T{1}) eslog::error("only in-place supported now (A=X, alpha=1)\n");

    size_t group_size = 64;
    sycl::nd_range<1> range(sycl::range<1>(B->get_size_primary() * group_size), sycl::range<1>(group_size));
    I * B_ptrs = B->ptrs;
    I * B_idxs = B->idxs;
    T * B_vals = B->vals;
    T * A_vals = A->vals;
    size_t A_ld = A->ld;
    T beta = this->beta;

    if(X->order == B->order) {
        q->q.parallel_for(range, [=](sycl::nd_item<1> item){
            sycl::group g = item.get_group();
            size_t ip = g.get_group_linear_id();
            T * A_ptr = A_vals + ip * A_ld;
            I start = B_ptrs[ip];
            I end = B_ptrs[ip+1];
            for(I i = start + g.get_local_linear_id(); i < end; i += g.get_local_linear_range()) {
                I is = B_idxs[i];
                T val = B_vals[i];
                A_ptr[is] += beta * val;
            }
        });
    }
    else {
        q->q.parallel_for(range, [=](sycl::nd_item<1> item){
            sycl::group g = item.get_group();
            size_t ip = g.get_group_linear_id();
            I start = B_ptrs[ip];
            I end = B_ptrs[ip+1];
            for(I i = start + g.get_local_linear_id(); i < end; i += g.get_local_linear_range()) {
                I is = B_idxs[i];
                T val = B_vals[i];
                A_vals[is * A_ld + ip] += beta * val;
            }
        });
    }
}



#define INSTANTIATE_T_I(T,I) \
template class w_oneapi_lincomb_ddnx_dcsy<T,I>;

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
