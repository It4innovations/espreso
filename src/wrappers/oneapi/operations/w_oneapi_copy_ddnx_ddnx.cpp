
#ifdef HAVE_ONEAPI

#include "wrappers/oneapi/operations/w_oneapi_copy_ddnx_ddnx.h"

#include "wrappers/oneapi/common_oneapi_mgm.h"
#include "wrappers/oneapi/common_oneblas.h"



namespace espreso {
namespace gpu {
namespace operations {



template<typename T>
void w_oneapi_copy_ddnx_ddnx<T>::internal_setup()
{
    wss_tmp_perform = 0;
}



template<typename T>
void w_oneapi_copy_ddnx_ddnx<T>::internal_perform(void * ws_tmp)
{
    if(uplo == 'L' || uplo == 'U') {
        size_t size_primary = M_src->get_size_primary();
        size_t size_secdary = M_src->get_size_secdary();
        T * src = M_src->vals;
        T * dst = M_dst->vals;
        size_t ld_src = M_src->ld;
        size_t ld_dst = M_dst->ld;
        int workitems_in_workgroup = 256;
        sycl::nd_range<1> range(sycl::range<1>(size_primary * workitems_in_workgroup), sycl::range<1>(workitems_in_workgroup));

        if((uplo == 'L') == (M_src->order == 'R')) {
            // copy_popullatedfirst
            q->q.parallel_for(range, [=](sycl::nd_item<1> item){
                sycl::group g = item.get_group();
                size_t ip = g.get_group_linear_id();
                T * sub_src = src + ld_src * ip;
                T * sub_dst = dst + ld_dst * ip;
                for(size_t is = g.get_local_linear_id(); is <= ip; is += g.get_local_linear_range()) {
                    sub_dst[is] = sub_src[is];
                }
            });
        }
        else {
            // copy_emptyfirst
            q->q.parallel_for(range, [=](sycl::nd_item<1> item){
                sycl::group g = item.get_group();
                size_t warp_size = item.get_sub_group().get_max_local_range().get(0);
                size_t ip = g.get_group_linear_id();
                T * sub_src = src + ld_src * ip;
                T * sub_dst = dst + ld_dst * ip;
                size_t start = (ip / warp_size) * warp_size;
                for(size_t is = start + g.get_local_linear_id(); is < size_secdary; is += g.get_local_linear_range()) {
                    if(is >= ip) {
                        sub_dst[is] = sub_src[is];
                    }
                }
            });
        }
    }
    else {
        oneblas::row_major::omatcopy(q->q, onemkl::transpose::nontrans, M_src->get_size_primary(), M_src->get_size_secdary(), T{1}, M_src->vals, M_src->ld, M_dst->vals, M_dst->ld);
    }
}



#define INSTANTIATE_T(T) \
template class w_oneapi_copy_ddnx_ddnx<T>;

    #define INSTANTIATE \
    /* INSTANTIATE_T(float) */ \
    INSTANTIATE_T(double) \
    /* INSTANTIATE_T(std::complex<float>) */ \
    INSTANTIATE_T(std::complex<double>)

        INSTANTIATE

    #undef INSTANTIATE
#undef INSTANTIATE_T



}
}
}

#endif
