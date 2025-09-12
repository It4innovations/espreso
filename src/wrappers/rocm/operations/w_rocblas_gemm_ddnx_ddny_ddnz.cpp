
#ifdef HAVE_ROCM

#include "wrappers/rocm/operations/w_rocblas_gemm_ddnx_ddny_ddnz.h"

#include "wrappers/rocm/common_rocblas.h"
#include "wrappers/rocm/common_rocm_mgm.h"



namespace espreso {
namespace gpu {
namespace operations {



struct w_rocblas_gemm_ddnx_ddny_ddnz_data
{
    rocblas_operation op_A;
    rocblas_operation op_B;
    bool swap_a_b;
    size_t m;
    size_t n;
    size_t k;
};



template<typename T>
w_rocblas_gemm_ddnx_ddny_ddnz<T>::w_rocblas_gemm_ddnx_ddny_ddnz() {}



template<typename T>
w_rocblas_gemm_ddnx_ddny_ddnz<T>::~w_rocblas_gemm_ddnx_ddny_ddnz() {}



template<typename T>
void w_rocblas_gemm_ddnx_ddny_ddnz<T>::internal_setup()
{
    data = std::make_unique<w_rocblas_gemm_ddnx_ddny_ddnz_data>();

    data->swap_a_b = (C->order == 'R');
    data->op_A = (((A->order == 'C') == (C->order == 'C')) ? rocblas_operation_none : rocblas_operation_transpose);
    data->op_B = (((B->order == 'C') == (C->order == 'C')) ? rocblas_operation_none : rocblas_operation_transpose);
    data->m = C->nrows;
    data->n = C->ncols;
    data->k = A->ncols;
    if(data->swap_a_b) {
        std::swap(data->m, data->n);
        std::swap(data->op_A, data->op_B);
    }

    CHECK(rocblas_start_device_memory_size_query(handle_dnblas->h));

    do_call('B');

    CHECK(rocblas_stop_device_memory_size_query(handle_dnblas->h, &wss_tmp_perform));
}



template<typename T>
void w_rocblas_gemm_ddnx_ddny_ddnz<T>::internal_perform(void * ws_tmp)
{
    CHECK(rocblas_set_workspace(handle_dnblas->h, ws_tmp, wss_tmp_perform));

    do_call('C');

    CHECK(rocblas_set_workspace(handle_dnblas->h, nullptr, 0));
}



template<typename T>
void w_rocblas_gemm_ddnx_ddny_ddnz<T>::do_call(char stage)
{
    using U = cpp_to_rocblas_type_t<T>;
    U * dummyptr_U = (U*)(sizeof(U));
    U * ptr_A_vals = ((stage == 'B') ? dummyptr_U : (U*)A->vals);
    U * ptr_B_vals = ((stage == 'B') ? dummyptr_U : (U*)B->vals);
    U * ptr_C_vals = ((stage == 'B') ? dummyptr_U : (U*)C->vals);
    size_t A_ld = A->ld;
    size_t B_ld = B->ld;
    size_t C_ld = C->ld;
    if(data->swap_a_b) {
        std::swap(ptr_A_vals, ptr_B_vals);
        std::swap(A_ld, B_ld);
    }

    if constexpr(std::is_same_v<T,float>)                CHECK(rocblas_sgemm(handle_dnblas->h, data->op_A, data->op_B, data->m, data->n, data->k, (U*)&alpha, ptr_A_vals, A_ld, ptr_B_vals, B_ld, (U*)&beta, ptr_C_vals, C_ld));
    if constexpr(std::is_same_v<T,double>)               CHECK(rocblas_dgemm(handle_dnblas->h, data->op_A, data->op_B, data->m, data->n, data->k, (U*)&alpha, ptr_A_vals, A_ld, ptr_B_vals, B_ld, (U*)&beta, ptr_C_vals, C_ld));
    if constexpr(std::is_same_v<T,std::complex<float>>)  CHECK(rocblas_cgemm(handle_dnblas->h, data->op_A, data->op_B, data->m, data->n, data->k, (U*)&alpha, ptr_A_vals, A_ld, ptr_B_vals, B_ld, (U*)&beta, ptr_C_vals, C_ld));
    if constexpr(std::is_same_v<T,std::complex<double>>) CHECK(rocblas_zgemm(handle_dnblas->h, data->op_A, data->op_B, data->m, data->n, data->k, (U*)&alpha, ptr_A_vals, A_ld, ptr_B_vals, B_ld, (U*)&beta, ptr_C_vals, C_ld));
}



#define INSTANTIATE_T(T) \
template class w_rocblas_gemm_ddnx_ddny_ddnz<T>;

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
