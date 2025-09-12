
#ifdef HAVE_ROCM

#include "wrappers/rocm/operations/w_rocblas_herk_ddnx_ddny.h"

#include "wrappers/rocm/common_rocblas.h"
#include "wrappers/rocm/common_rocm_mgm.h"



namespace espreso {
namespace gpu {
namespace operations {



struct w_rocblas_herk_ddnx_ddny_data
{
    rocblas_fill uplo_C;
    rocblas_operation trans_mode;
    size_t n;
    size_t k;
};



template<typename T>
w_rocblas_herk_ddnx_ddny<T>::w_rocblas_herk_ddnx_ddny() {}



template<typename T>
w_rocblas_herk_ddnx_ddny<T>::~w_rocblas_herk_ddnx_ddny() {}



template<typename T>
void w_rocblas_herk_ddnx_ddny<T>::internal_setup()
{
    data = std::make_unique<w_rocblas_herk_ddnx_ddny_data>();

    data->uplo_C = (((C->prop.uplo == 'L') == (C->order == 'C')) ? rocblas_fill_lower : rocblas_fill_upper);
    data->trans_mode = (((mode == math::blas::herk_mode::AAh) == (A->order == 'C')) ? rocblas_operation_none : rocblas_operation_transpose);
    data->n = A->nrows;
    data->k = A->ncols;
    if(mode == math::blas::herk_mode::AhA) {
        std::swap(data->n, data->k);
    }

    CHECK(rocblas_start_device_memory_size_query(handle_dnblas->h));

    do_call('B');

    CHECK(rocblas_stop_device_memory_size_query(handle_dnblas->h, &wss_tmp_perform));
}



template<typename T>
void w_rocblas_herk_ddnx_ddny<T>::internal_perform(void * ws_tmp)
{
    CHECK(rocblas_set_workspace(handle_dnblas->h, ws_tmp, wss_tmp_perform));

    do_call('C');

    CHECK(rocblas_set_workspace(handle_dnblas->h, nullptr, 0));
}



template<typename T>
void w_rocblas_herk_ddnx_ddny<T>::do_call(char stage)
{
    if(utils::is_complex<T>()) eslog::error("complex types not yet supported in herk\n");
    using U = cpp_to_rocblas_type_t<T>;
    U * dummyptr_U = (U*)(sizeof(U));
    U * ptr_A_vals = ((stage == 'B') ? dummyptr_U : (U*)A->vals);
    U * ptr_C_vals = ((stage == 'B') ? dummyptr_U : (U*)C->vals);

    if constexpr(std::is_same_v<T,float>)                CHECK(rocblas_ssyrk(handle_dnblas->h, data->uplo_C, data->trans_mode, data->n, data->k, &alpha, ptr_A_vals, A->ld, &beta, ptr_C_vals, C->ld));
    if constexpr(std::is_same_v<T,double>)               CHECK(rocblas_dsyrk(handle_dnblas->h, data->uplo_C, data->trans_mode, data->n, data->k, &alpha, ptr_A_vals, A->ld, &beta, ptr_C_vals, C->ld));
    if constexpr(std::is_same_v<T,std::complex<float>>)  CHECK(rocblas_cherk(handle_dnblas->h, data->uplo_C, data->trans_mode, data->n, data->k, &alpha, ptr_A_vals, A->ld, &beta, ptr_C_vals, C->ld));
    if constexpr(std::is_same_v<T,std::complex<double>>) CHECK(rocblas_zherk(handle_dnblas->h, data->uplo_C, data->trans_mode, data->n, data->k, &alpha, ptr_A_vals, A->ld, &beta, ptr_C_vals, C->ld));
}



#define INSTANTIATE_T(T) \
template class w_rocblas_herk_ddnx_ddny<T>;

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
