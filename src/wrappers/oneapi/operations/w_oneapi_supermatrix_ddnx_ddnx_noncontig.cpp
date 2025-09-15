
#ifdef HAVE_ONEAPI

#include "wrappers/oneapi/operations/w_oneapi_supermatrix_ddnx_ddnx_noncontig.h"

#include "wrappers/oneapi/common_oneapi_mgm.h"
#include "wrappers/oneapi/common_oneblas.h"
#include "wrappers/oneapi/common_internal.hpp"



namespace espreso {
namespace gpu {
namespace operations {



template<typename T, typename I, bool SUP_PRIMARY, bool SUP_SECDARY, typename OP>
struct functor_do_supermatrix
{
    I size_src_primary;
    I size_src_secdary;
    T * src;
    I ld_src;
    T * dst;
    I ld_dst;
    I * map_primary;
    I * map_secdary;
    OP op;
    functor_do_supermatrix(I size_src_primary_, I size_src_secdary_, T * src_, I ld_src_, T * dst_, I ld_dst_, I * map_primary_, I * map_secdary_, OP op_)
        : size_src_primary(size_src_primary_)
        , size_src_secdary(size_src_secdary_)
        , src(src_)
        , ld_src(ld_src_)
        , dst(dst_)
        , ld_dst(ld_dst_)
        , map_primary(map_primary_)
        , map_secdary(map_secdary_)
        , op(op_)
    {}
    void operator()(sycl::nd_item<2> item) const {
        sycl::group g = item.get_group();
        I ips_start = g.get_group_id(0) * g.get_local_range(0) + g.get_local_id(0);
        I ips_stride = g.get_local_range(0) * g.get_group_range(0);
        I iss_start = g.get_group_id(1) * g.get_local_range(1) + g.get_local_id(1);
        I iss_stride = g.get_local_range(1) * g.get_group_range(1);

        for(I ips = ips_start; ips < size_src_primary; ips += ips_stride) {
            I ipd = ips;
            if constexpr(SUP_PRIMARY) {
                ipd = map_primary[ips];
            }

            T * dst_prim = dst + ipd * ld_dst;
            T * src_prim = src + ips * ld_src;

            for(I iss = iss_start; iss < size_src_secdary; iss += iss_stride) {
                I isd = iss;
                if constexpr(SUP_SECDARY) {
                    isd = map_secdary[iss];
                }

                op(dst_prim[isd], src_prim[iss]);
            }
        }
    }
};



template<typename T, typename I, typename OP>
static void do_supermatrix_launch_kernel_2(sycl::queue & q, size_t size_src_primary, size_t size_src_secdary, T * src, size_t ld_src, T * dst, size_t ld_dst, I * map_primary, I * map_secdary, OP op)
{
    sycl::range local_range(8,32);
    sycl::range n(2,2);
    sycl::range global_range(utils::divide_round_up(size_src_primary, n[0]), utils::divide_round_up(size_src_secdary, n[1]));
    sycl::nd_range range(global_range, local_range);

    if(map_primary == nullptr && map_secdary == nullptr) {
        q.parallel_for(range, functor_do_supermatrix<T,I,false,false,OP>(size_src_primary, size_src_secdary, src, ld_src, dst, ld_dst, map_primary, map_secdary, op));
    }
    if(map_primary != nullptr && map_secdary == nullptr) {
        q.parallel_for(range, functor_do_supermatrix<T,I,true,false,OP>(size_src_primary, size_src_secdary, src, ld_src, dst, ld_dst, map_primary, map_secdary, op));
    }
    if(map_primary == nullptr && map_secdary != nullptr) {
        q.parallel_for(range, functor_do_supermatrix<T,I,false,true,OP>(size_src_primary, size_src_secdary, src, ld_src, dst, ld_dst, map_primary, map_secdary, op));
    }
    if(map_primary != nullptr && map_secdary != nullptr) {
        q.parallel_for(range, functor_do_supermatrix<T,I,true,true,OP>(size_src_primary, size_src_secdary, src, ld_src, dst, ld_dst, map_primary, map_secdary, op));
    }
}



template<typename T, typename I>
static void do_supermatrix_launch_kernel_1(sycl::queue & q, size_t size_src_primary, size_t size_src_secdary, T * src, size_t ld_src, T * dst, size_t ld_dst, I * map_primary, I * map_secdary, typename supermatrix_ddnx_ddnx_noncontig<T,I>::mode mode_val)
{
    using mode = typename supermatrix_ddnx_ddnx_noncontig<T,I>::mode;

    switch(mode_val) {
        case mode::assign: {
            auto op = [] (T & dst, T src){ dst = src; };
            do_supermatrix_launch_kernel_2(q, size_src_primary, size_src_secdary, src, ld_src, dst, ld_dst, map_primary, map_secdary, op);
            break;
        }
        case mode::accumulate: {
            auto op = [] (T & dst, T src){ dst += src; };
            do_supermatrix_launch_kernel_2(q, size_src_primary, size_src_secdary, src, ld_src, dst, ld_dst, map_primary, map_secdary, op);
            break;
        }
        case mode::accumulate_atomic: {
            auto op = [] (T & dst, T src){ my_atomicadd(&dst, src); };
            do_supermatrix_launch_kernel_2(q, size_src_primary, size_src_secdary, src, ld_src, dst, ld_dst, map_primary, map_secdary, op);
            break;
        }
        default: {
            eslog::error("wrong mode\n");
        }
    }
}



template<typename T, typename I>
void w_oneapi_supermatrix_ddnx_ddnx_noncontig<T,I>::internal_perform()
{
    VectorDenseView_new<I> * map_primary = ((d_M_src->order == 'R') ? d_row_map : d_col_map);
    VectorDenseView_new<I> * map_secdary = ((d_M_src->order == 'R') ? d_col_map : d_row_map);
    I * map_primary_vals = ((map_primary != nullptr) ? map_primary->vals : nullptr);
    I * map_secdary_vals = ((map_secdary != nullptr) ? map_secdary->vals : nullptr);

    do_supermatrix_launch_kernel_1<T,I>(q->q, d_M_src->get_size_primary(), d_M_src->get_size_secdary(), d_M_src->vals, d_M_src->ld, d_M_dst->vals, d_M_dst->ld, map_primary_vals, map_secdary_vals, mode_val);
}



#define INSTANTIATE_T_I(T,I) \
template class w_oneapi_supermatrix_ddnx_ddnx_noncontig<T,I>;

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
