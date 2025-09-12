
#ifdef HAVE_ROCM

#include "wrappers/rocm/operations/w_rocsparse_trsm_dcsx_ddny_ddny.h"

#include "wrappers/rocm/common_rocsparse.h"
#include "wrappers/rocm/common_rocm_mgm.h"



namespace espreso {
namespace gpu {
namespace operations {



struct w_rocsparse_trsm_dcsx_ddny_ddny_data
{
    rocsparse_spmat_descr descr_A;
    rocsparse_dnmat_descr descr_X;
    rocsparse_dnmat_descr descr_B;
    rocsparse_operation op_A;
    rocsparse_operation op_B;
    rocsparse_spsm_alg alg;
};



template<typename T, typename I>
w_rocsparse_trsm_dcsx_ddny_ddny<T,I>::w_rocsparse_trsm_dcsx_ddny_ddny() {}



template<typename T, typename I>
w_rocsparse_trsm_dcsx_ddny_ddny<T,I>::~w_rocsparse_trsm_dcsx_ddny_ddny()
{
    if(this->called_setup) {
        CHECK(rocsparse_destroy_spmat_descr(data->descr_A));
        CHECK(rocsparse_destroy_dnmat_descr(data->descr_X));
        CHECK(rocsparse_destroy_dnmat_descr(data->descr_B));
    }
}



template<typename T, typename I>
char w_rocsparse_trsm_dcsx_ddny_ddny<T,I>::internal_get_native_place()
{
    return 'O';
}



template<typename T, typename I>
void w_rocsparse_trsm_dcsx_ddny_ddny<T,I>::internal_setup()
{
    data = std::make_unique<w_rocsparse_trsm_dcsx_ddny_ddny_data>();

#if HIP_VERSION_MAJOR > 6 || (HIP_VERSION_MAJOR == 6 && HIP_VERSION_MINOR >= 2)
    static_assert(false, "not sure how (not) buggy it is in newer versions, check and redo");
#else
    // older rocsparse assumes colmajor for both dense matrices
    // older rocsparse mistakenly thinks that trans_B is applied to both B and C
    // internally, it is performed in-place. B is copied into C, then in-place with C

    if(B->order != X->order) eslog::error("in rocsparse, B and X must have equal order\n");

    T * dummyptr_T = (T*)(sizeof(T));
    I * dummyptr_I = (I*)(sizeof(I));
    if(A->order == 'R') CHECK(rocsparse_create_csr_descr(&data->descr_A, A->nrows, A->ncols, A->nnz, dummyptr_I, dummyptr_I, dummyptr_T, get_rocsparse_index_type<I>(), get_rocsparse_index_type<I>(), rocsparse_index_base_zero, get_rocsparse_data_type<T>()));
    if(A->order == 'C') CHECK(rocsparse_create_csr_descr(&data->descr_A, A->ncols, A->nrows, A->nnz, dummyptr_I, dummyptr_I, dummyptr_T, get_rocsparse_index_type<I>(), get_rocsparse_index_type<I>(), rocsparse_index_base_zero, get_rocsparse_data_type<T>()));
    if(B->order == 'R') CHECK(rocsparse_create_dnmat_descr(&data->descr_B, B->ncols, B->nrows, B->ld, dummyptr_T, get_rocsparse_data_type<T>(), rocsparse_order_column));
    if(B->order == 'C') CHECK(rocsparse_create_dnmat_descr(&data->descr_B, B->nrows, B->ncols, B->ld, dummyptr_T, get_rocsparse_data_type<T>(), rocsparse_order_column));
    if(X->order == 'R') CHECK(rocsparse_create_dnmat_descr(&data->descr_X, X->ncols, X->nrows, X->ld, dummyptr_T, get_rocsparse_data_type<T>(), rocsparse_order_column));
    if(X->order == 'C') CHECK(rocsparse_create_dnmat_descr(&data->descr_X, X->nrows, X->ncols, X->ld, dummyptr_T, get_rocsparse_data_type<T>(), rocsparse_order_column));

    rocsparse_fill_mode upper = rocsparse_fill_mode_upper;
    rocsparse_fill_mode lower = rocsparse_fill_mode_lower;
    rocsparse_diag_type unit = rocsparse_diag_type_unit;
    rocsparse_diag_type nonunit = rocsparse_diag_type_non_unit;
    // rocsparse_matrix_type triangular = rocsparse_matrix_type_triangular;
    rocsparse_storage_mode sorted = rocsparse_storage_mode_sorted;
    if(A->prop.uplo == 'L') CHECK(rocsparse_spmat_set_attribute(data->descr_A, rocsparse_spmat_fill_mode, &lower, sizeof(lower)));
    if(A->prop.uplo == 'U') CHECK(rocsparse_spmat_set_attribute(data->descr_A, rocsparse_spmat_fill_mode, &upper, sizeof(upper)));
    if(A->prop.diag == 'U') CHECK(rocsparse_spmat_set_attribute(data->descr_A, rocsparse_spmat_diag_type, &unit, sizeof(unit)));
    if(A->prop.diag == 'N') CHECK(rocsparse_spmat_set_attribute(data->descr_A, rocsparse_spmat_diag_type, &nonunit, sizeof(nonunit)));
    // if(fill != 'N') CHECK(rocsparse_spmat_set_attribute(data->descr_A, rocsparse_spmat_matrix_type, &triangular, sizeof(triangular)));
    CHECK(rocsparse_spmat_set_attribute(data->descr_A, rocsparse_spmat_storage_mode, &sorted, sizeof(sorted)));

    data->op_A = ((A->order == 'R') ? rocsparse_operation_none : rocsparse_operation_transpose);
    data->op_B = ((X->order == 'C') ? rocsparse_operation_none : rocsparse_operation_transpose);

    data->alg = rocsparse_spsm_alg_default;

    size_t buffersize;
    T alpha = T{1};
    CHECK(rocsparse_spsm(handle_spblas->h, data->op_A, data->op_B, &alpha, data->descr_A, data->descr_B, data->descr_X, get_rocsparse_data_type<T>(), data->alg, rocsparse_spsm_stage_buffer_size, &buffersize, nullptr));

    wss_internal = 0;
    wss_persistent = 0;
    wss_tmp_preprocess = buffersize;
    wss_tmp_perform = buffersize;
#endif
}



template<typename T, typename I>
void w_rocsparse_trsm_dcsx_ddny_ddny<T,I>::internal_preprocess(void * ws_tmp)
{
#if HIP_VERSION_MAJOR > 6 || (HIP_VERSION_MAJOR == 6 && HIP_VERSION_MINOR >= 2)
    static_assert(false, "not sure how (not) buggy it is in newer versions, check and redo");
#else
    if(A->order == 'R') CHECK(rocsparse_csr_set_pointers(data->descr_A, A->ptrs, A->idxs, A->vals));
    if(A->order == 'C') CHECK(rocsparse_csc_set_pointers(data->descr_A, A->ptrs, A->idxs, A->vals));

    T alpha = T{1};
    CHECK(rocsparse_spsm(handle_spblas->h, data->op_A, data->op_B, &alpha, data->descr_A, data->descr_B, data->descr_X, get_rocsparse_data_type<T>(), data->alg, rocsparse_spsm_stage_preprocess, &wss_tmp_preprocess, ws_tmp));
#endif


}



template<typename T, typename I>
void w_rocsparse_trsm_dcsx_ddny_ddny<T,I>::internal_perform(void * ws_tmp)
{
#if HIP_VERSION_MAJOR > 6 || (HIP_VERSION_MAJOR == 6 && HIP_VERSION_MINOR >= 2)
    static_assert(false, "not sure how (not) buggy it is in newer versions, check and redo");
#else
    if(A->order == 'R') CHECK(rocsparse_csr_set_pointers(data->descr_A, A->ptrs, A->idxs, A->vals));
    if(A->order == 'C') CHECK(rocsparse_csc_set_pointers(data->descr_A, A->ptrs, A->idxs, A->vals));
    CHECK(rocsparse_dnmat_set_values(data->descr_B, B->vals));
    CHECK(rocsparse_dnmat_set_values(data->descr_X, X->vals));

    // no update of internal matrix structure numerical values needed

    T alpha = T{1};
    CHECK(rocsparse_spsm(handle_spblas->h, data->op_A, data->op_B, &alpha, data->descr_A, data->descr_B, data->descr_X, get_rocsparse_data_type<T>(), data->alg, rocsparse_spsm_stage_compute, &wss_tmp_perform, ws_tmp));

#endif


}



#define INSTANTIATE_T_I(T,I) \
template class w_rocsparse_trsm_dcsx_ddny_ddny<T,I>;

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
