
#ifdef HAVE_ROCM

#include "wrappers/rocm/operations/w_rocblas_trsm_ddnx_ddny.h"

#include "wrappers/rocm/common_rocblas.h"



namespace espreso {
namespace gpu {
namespace operations {



struct w_rocblas_trsm_ddnx_ddny_data
{
    rocblas_side side;
    rocblas_fill uplo_A;
    rocblas_operation op_A;
    rocblas_diagonal diag_A;
    size_t m;
    size_t n;
};



template<typename T>
w_rocblas_trsm_ddnx_ddny<T>::w_rocblas_trsm_ddnx_ddny() {}



template<typename T>
w_rocblas_trsm_ddnx_ddny<T>::~w_rocblas_trsm_ddnx_ddny() {}



template<typename T>
void w_rocblas_trsm_ddnx_ddny<T>::internal_setup()
{
    data = std::make_unique<w_rocblas_trsm_ddnx_ddny_data>();

    data->side = ((X->order == 'C') ? rocblas_side_left : rocblas_side_right);
    data->op_A = ((X->order == A->order) ? rocblas_operation_none : rocblas_operation_transpose);
    data->uplo_A = (((A->prop.uplo == 'L') == (A->order == 'C')) ? rocblas_fill_lower : rocblas_fill_upper);
    data->diag_A = ((A->prop.diag == 'U') ? rocblas_diagonal_unit : rocblas_diagonal_non_unit);
    data->m = X->nrows;
    data->n = X->ncols;
    if(X->order == 'R') {
        std::swap(data->m, data->n);
    }

    CHECK(rocblas_start_device_memory_size_query(handle_dnblas->h));

    do_call('B');

    CHECK(rocblas_stop_device_memory_size_query(handle_dnblas->h, &wss_tmp_perform));
}



template<typename T>
void w_rocblas_trsm_ddnx_ddny<T>::internal_perform(void * ws_tmp)
{
    CHECK(rocblas_set_workspace(handle_dnblas->h, ws_tmp, wss_tmp_perform));

    do_call('C');

    CHECK(rocblas_set_workspace(handle_dnblas->h, nullptr, 0));
}



template<typename T>
void w_rocblas_trsm_ddnx_ddny<T>::do_call(char stage)
{
    using U = cpp_to_rocblas_type_t<T>;
    U one = U{1};
    U * dummyptr_U = (U*)(sizeof(U));
    U * ptr_A_vals = ((stage == 'B') ? dummyptr_U : (U*)A->vals);
    U * ptr_X_vals = ((stage == 'B') ? dummyptr_U : (U*)X->vals);
    
    if constexpr(std::is_same_v<T,float>)                CHECK(rocblas_strsm(handle_dnblas->h, data->side, data->uplo_A, data->op_A, data->diag_A, data->m, data->n, &one, ptr_A_vals, A->ld, ptr_X_vals, X->ld));
    if constexpr(std::is_same_v<T,double>)               CHECK(rocblas_dtrsm(handle_dnblas->h, data->side, data->uplo_A, data->op_A, data->diag_A, data->m, data->n, &one, ptr_A_vals, A->ld, ptr_X_vals, X->ld));
    if constexpr(std::is_same_v<T,std::complex<float>>)  CHECK(rocblas_ctrsm(handle_dnblas->h, data->side, data->uplo_A, data->op_A, data->diag_A, data->m, data->n, &one, ptr_A_vals, A->ld, ptr_X_vals, X->ld));
    if constexpr(std::is_same_v<T,std::complex<double>>) CHECK(rocblas_ztrsm(handle_dnblas->h, data->side, data->uplo_A, data->op_A, data->diag_A, data->m, data->n, &one, ptr_A_vals, A->ld, ptr_X_vals, X->ld));
}



#define INSTANTIATE_T(T) \
template class w_rocblas_trsm_ddnx_ddny<T>;

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
