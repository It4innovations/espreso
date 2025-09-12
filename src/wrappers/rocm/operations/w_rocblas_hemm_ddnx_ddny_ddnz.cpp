
#ifdef HAVE_ROCM

#include "wrappers/rocm/operations/w_rocblas_hemm_ddnx_ddny_ddnz.h"

#include "wrappers/rocm/common_rocblas.h"
#include "wrappers/rocm/common_rocm_mgm.h"



namespace espreso {
namespace gpu {
namespace operations {



struct w_rocblas_hemm_ddnx_ddny_ddnz_data
{
    rocblas_side side;
    rocblas_fill uplo;
    size_t m;
    size_t n;
};



template<typename T>
w_rocblas_hemm_ddnx_ddny_ddnz<T>::w_rocblas_hemm_ddnx_ddny_ddnz() {}



template<typename T>
w_rocblas_hemm_ddnx_ddny_ddnz<T>::~w_rocblas_hemm_ddnx_ddny_ddnz() {}



template<typename T>
void w_rocblas_hemm_ddnx_ddny_ddnz<T>::internal_setup()
{
    data = std::make_unique<w_rocblas_hemm_ddnx_ddny_ddnz_data>();

    if(B->order != C->order) eslog::error("order of matrices B and C must match\n");
    if(utils::is_complex<T>() && A->order != C->order) eslog::error("for complex, order of all matrices must match\n");

    data->side = ((C->order == 'C') ? rocblas_side_left : rocblas_side_right);
    data->uplo = ((((C->order == 'C') == (A->order == 'C')) == (A->prop.uplo == 'L')) ? rocblas_fill_lower : rocblas_fill_upper);

    data->m = C->nrows;
    data->n = C->ncols;
    if(data->side == rocblas_side_right) std::swap(data->m, data->n);

    CHECK(rocblas_start_device_memory_size_query(handle_dnblas->h));

    do_call('B');

    CHECK(rocblas_stop_device_memory_size_query(handle_dnblas->h, &wss_tmp_perform));
}



template<typename T>
void w_rocblas_hemm_ddnx_ddny_ddnz<T>::internal_perform(void * ws_tmp)
{
    CHECK(rocblas_set_workspace(handle_dnblas->h, ws_tmp, wss_tmp_perform));

    do_call('C');

    CHECK(rocblas_set_workspace(handle_dnblas->h, nullptr, 0));
}



template<typename T>
void w_rocblas_hemm_ddnx_ddny_ddnz<T>::do_call(char stage)
{
    using U = cpp_to_rocblas_type_t<T>;
    U * dummyptr_U = (U*)(sizeof(U));
    U * ptr_A_vals = ((stage == 'B') ? dummyptr_U : (U*)A->vals);
    U * ptr_B_vals = ((stage == 'B') ? dummyptr_U : (U*)B->vals);
    U * ptr_C_vals = ((stage == 'B') ? dummyptr_U : (U*)C->vals);

    if constexpr(std::is_same_v<T,float>)                CHECK(rocblas_ssymm(handle_dnblas->h, data->side, data->uplo, data->m, data->n, (U*)&alpha, ptr_A_vals, A->ld, ptr_B_vals, B->ld, (U*)&beta, ptr_C_vals, C->ld));
    if constexpr(std::is_same_v<T,double>)               CHECK(rocblas_dsymm(handle_dnblas->h, data->side, data->uplo, data->m, data->n, (U*)&alpha, ptr_A_vals, A->ld, ptr_B_vals, B->ld, (U*)&beta, ptr_C_vals, C->ld));
    if constexpr(std::is_same_v<T,std::complex<float>>)  CHECK(rocblas_chemm(handle_dnblas->h, data->side, data->uplo, data->m, data->n, (U*)&alpha, ptr_A_vals, A->ld, ptr_B_vals, B->ld, (U*)&beta, ptr_C_vals, C->ld));
    if constexpr(std::is_same_v<T,std::complex<double>>) CHECK(rocblas_zhemm(handle_dnblas->h, data->side, data->uplo, data->m, data->n, (U*)&alpha, ptr_A_vals, A->ld, ptr_B_vals, B->ld, (U*)&beta, ptr_C_vals, C->ld));
}



#define INSTANTIATE_T(T) \
template class w_rocblas_hemm_ddnx_ddny_ddnz<T>;

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
