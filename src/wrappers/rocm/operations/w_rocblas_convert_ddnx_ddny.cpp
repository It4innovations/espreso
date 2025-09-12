
#ifdef HAVE_ROCM

#include "wrappers/rocm/operations/w_rocblas_convert_ddnx_ddny.h"

#include "wrappers/rocm/common_rocm_mgm.h"
#include "wrappers/rocm/common_rocblas.h"



namespace espreso {
namespace gpu {
namespace operations {



template<typename T>
void w_rocblas_convert_ddnx_ddny<T>::internal_setup()
{
    CHECK(rocblas_start_device_memory_size_query(handle_dnblas->h));

    do_call('B');

    CHECK(rocblas_stop_device_memory_size_query(handle_dnblas->h, &wss_tmp_perform));
}



template<typename T>
void w_rocblas_convert_ddnx_ddny<T>::internal_perform(void * ws_tmp)
{
    CHECK(rocblas_set_workspace(handle_dnblas->h, ws_tmp, wss_tmp_perform));

    do_call('C');

    CHECK(rocblas_set_workspace(handle_dnblas->h, nullptr, 0));
}



template<typename T>
void w_rocblas_convert_ddnx_ddny<T>::do_call(char stage)
{
    if(M_src->order == M_dst->order) {
        if(stage == 'C') {
            CHECK(hipMemcpy2DAsync(M_dst->vals, M_dst->ld * sizeof(T), M_src->vals, M_src->ld * sizeof(T), M_src->get_size_secdary() * sizeof(T), M_src->get_size_primary(), hipMemcpyDefault, q->stream));
        }
    }
    else {
        using U = cpp_to_rocblas_type_t<T>;
        T one = T{1};
        T zero = T{0};
        U * dummyptr_U = (U*)(sizeof(U));
        U * ptr_src_vals = ((stage == 'B') ? dummyptr_U : (U*)M_src->vals);
        U * ptr_dst_vals = ((stage == 'B') ? dummyptr_U : (U*)M_dst->vals);

        if constexpr(std::is_same_v<T,float>)                CHECK(rocblas_sgeam(handle_dnblas->h, rocblas_operation_transpose, rocblas_operation_none, M_src->get_size_secdary(), M_src->get_size_primary(), (U*)&one, ptr_src_vals, M_src->ld, (U*)&zero, nullptr, M_src->ld, ptr_dst_vals, M_dst->ld));
        if constexpr(std::is_same_v<T,double>)               CHECK(rocblas_dgeam(handle_dnblas->h, rocblas_operation_transpose, rocblas_operation_none, M_src->get_size_secdary(), M_src->get_size_primary(), (U*)&one, ptr_src_vals, M_src->ld, (U*)&zero, nullptr, M_src->ld, ptr_dst_vals, M_dst->ld));
        if constexpr(std::is_same_v<T,std::complex<float>>)  CHECK(rocblas_cgeam(handle_dnblas->h, rocblas_operation_transpose, rocblas_operation_none, M_src->get_size_secdary(), M_src->get_size_primary(), (U*)&one, ptr_src_vals, M_src->ld, (U*)&zero, nullptr, M_src->ld, ptr_dst_vals, M_dst->ld));
        if constexpr(std::is_same_v<T,std::complex<double>>) CHECK(rocblas_zgeam(handle_dnblas->h, rocblas_operation_transpose, rocblas_operation_none, M_src->get_size_secdary(), M_src->get_size_primary(), (U*)&one, ptr_src_vals, M_src->ld, (U*)&zero, nullptr, M_src->ld, ptr_dst_vals, M_dst->ld));
    }
}



#define INSTANTIATE_T(T) \
template class w_rocblas_convert_ddnx_ddny<T>;

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
