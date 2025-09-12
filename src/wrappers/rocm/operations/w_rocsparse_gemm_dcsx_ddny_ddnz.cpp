
#ifdef HAVE_ROCM

#include "wrappers/rocm/operations/w_rocsparse_gemm_dcsx_ddny_ddnz.h"

#include "wrappers/rocm/common_rocsparse.h"



namespace espreso {
namespace gpu {
namespace operations {



template<typename T, typename I>
struct w_rocsparse_gemm_dcsx_ddny_ddnz_data
{
    rocsparse_spmat_descr descr_A;
    rocsparse_dnmat_descr descr_B;
    rocsparse_dnmat_descr descr_C;
    rocsparse_spmm_alg alg;
};



template<typename T, typename I>
w_rocsparse_gemm_dcsx_ddny_ddnz<T,I>::w_rocsparse_gemm_dcsx_ddny_ddnz() {}



template<typename T, typename I>
w_rocsparse_gemm_dcsx_ddny_ddnz<T,I>::~w_rocsparse_gemm_dcsx_ddny_ddnz()
{
    if(this->called_setup) {
        CHECK(rocsparse_destroy_spmat_descr(data->descr_A));
        CHECK(rocsparse_destroy_dnmat_descr(data->descr_B));
        CHECK(rocsparse_destroy_dnmat_descr(data->descr_C));
    }
}



template<typename T, typename I>
void w_rocsparse_gemm_dcsx_ddny_ddnz<T,I>::internal_setup()
{
    data = std::make_unique<w_rocsparse_gemm_dcsx_ddny_ddnz_data<T,I>>();

    T * dummyptr_T = (T*)(sizeof(T));
    I * dummyptr_I = (I*)(sizeof(I));
    if(A->order == 'R') CHECK(rocsparse_create_csr_descr(&data->descr_A, A->nrows, A->ncols, A->nnz, dummyptr_I, dummyptr_I, dummyptr_T, get_rocsparse_index_type<I>(), get_rocsparse_index_type<I>(), rocsparse_index_base_zero, get_rocsparse_data_type<T>()));
    if(A->order == 'C') CHECK(rocsparse_create_csc_descr(&data->descr_A, A->nrows, A->ncols, A->nnz, dummyptr_I, dummyptr_I, dummyptr_T, get_rocsparse_index_type<I>(), get_rocsparse_index_type<I>(), rocsparse_index_base_zero, get_rocsparse_data_type<T>()));

    CHECK(rocsparse_create_dnmat_descr(&data->descr_B, B->nrows, B->ncols, B->ld, dummyptr_T, get_rocsparse_data_type<T>(), get_rocsparse_order(B->order)));
    CHECK(rocsparse_create_dnmat_descr(&data->descr_C, C->nrows, C->ncols, C->ld, dummyptr_T, get_rocsparse_data_type<T>(), get_rocsparse_order(C->order)));

    data->alg = rocsparse_spmm_alg_default;

    wss_internal = 0;
    CHECK(rocsparse_spmm(handle_spblas->h, rocsparse_operation_none, rocsparse_operation_none, &alpha, data->descr_A, data->descr_B, &beta, data->descr_C, get_rocsparse_data_type<T>(), data->alg, rocsparse_spmm_stage_buffer_size, &wss_persistent, nullptr));
    wss_tmp_preprocess = 0;
    wss_tmp_perform = 0;
}



template<typename T, typename I>
void w_rocsparse_gemm_dcsx_ddny_ddnz<T,I>::internal_preprocess(void * /*ws_tmp*/)
{
    if(A->order == 'R') CHECK(rocsparse_csr_set_pointers(data->descr_A, A->ptrs, A->idxs, A->vals));
    if(A->order == 'C') CHECK(rocsparse_csc_set_pointers(data->descr_A, A->ptrs, A->idxs, A->vals));

    CHECK(rocsparse_spmm(handle_spblas->h, rocsparse_operation_none, rocsparse_operation_none, &alpha, data->descr_A, data->descr_B, &beta, data->descr_C, get_rocsparse_data_type<T>(), data->alg, rocsparse_spmm_stage_preprocess, &wss_persistent, ws_persistent));
}



template<typename T, typename I>
void w_rocsparse_gemm_dcsx_ddny_ddnz<T,I>::internal_perform(void * /*ws_tmp*/)
{
    if(A->order == 'R') CHECK(rocsparse_csr_set_pointers(data->descr_A, A->ptrs, A->idxs, A->vals));
    if(A->order == 'C') CHECK(rocsparse_csc_set_pointers(data->descr_A, A->ptrs, A->idxs, A->vals));
    CHECK(rocsparse_dnmat_set_values(data->descr_B, B->vals));
    CHECK(rocsparse_dnmat_set_values(data->descr_C, C->vals));

    CHECK(rocsparse_spmm(handle_spblas->h, rocsparse_operation_none, rocsparse_operation_none, &alpha, data->descr_A, data->descr_B, &beta, data->descr_C, get_rocsparse_data_type<T>(), data->alg, rocsparse_spmm_stage_compute, &wss_persistent, ws_persistent));
}



#define INSTANTIATE_T_I(T,I) \
template class w_rocsparse_gemm_dcsx_ddny_ddnz<T,I>;

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
