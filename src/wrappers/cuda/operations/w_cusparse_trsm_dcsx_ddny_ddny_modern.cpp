
#ifdef HAVE_CUDA
#ifndef USE_CUSPARSE_LEGACY

#include "wrappers/cuda/operations/trsm_dcsx_ddny_ddny.h"

#include <cusparse.h>

#include "wrappers/cuda/common_cusparse.h"



namespace espreso {
namespace gpu {
namespace operations {



struct w_cusparse_trsm_dcsx_ddny_ddny_data
{
    cusparseHandle_t handle_cusparse;
    cusparseSpMatDescr_t descr_A;
    cusparseDnMatDescr_t descr_X;
    cusparseDnMatDescr_t descr_B;
    cusparseSpSMDescr_t descr_spsm;
    cusparseOperation_t op_A;
    cusparseOperation_t op_B;
    cusparseSpSMAlg_t spsm_alg;
};



template<typename T, typename I>
w_cusparse_trsm_dcsx_ddny_ddny<T,I>::w_cusparse_trsm_dcsx_ddny_ddny()
{
    data = std::make_unique<w_cusparse_trsm_dcsx_ddny_ddny_data>();
}



template<typename T, typename I>
w_cusparse_trsm_dcsx_ddny_ddny<T,I>::~w_cusparse_trsm_dcsx_ddny_ddny()
{
    if(called_set_A) {
        CHECK(cusparseDestroySpMat(data->descr_A));
    }
    if(called_set_X) {
        CHECK(cusparseDestroyDnMat(data->descr_X));
    }
    if(called_set_B) {
        CHECK(cusparseDestroyDnMat(data->descr_B));
    }
    if(called_setup) {
        CHECK(cusparseSpSM_destroyDescr(data->descr_spsm));
    }

    data.reset();
}



template<typename T, typename I>
char w_cusparse_trsm_dcsx_ddny_ddny<T,I>::internal_get_native_place()
{
    return 'I';
}



template<typename T, typename I>
void w_cusparse_trsm_dcsx_ddny_ddny<T,I>::internal_set_matrix_A()
{
    if(called_set_A) {
        CHECK(cusparseCsrSetPointers(data->descr_A, A.ptrs, A.idxs, A.vals));
    }
    else {
        CHECK(cusparseCreateCsr(data->descr_A, A.nrows, A.ncols, A.nnz, A.ptrs, A.idxs, A.vals, cusparse_index_type<I>(), cusparse_index_type<I>(), CUSPARSE_INDEX_BASE_ZERO, cusparse_data_type<T>()));
        auto nonunit = CUSPARSE_DIAG_TYPE_NON_UNIT;
        auto lower = CUSPARSE_FILL_MODE_LOWER;
        auto upper = CUSPARSE_FILL_MODE_UPPER;
        CHECK(cusparseSpMatSetAttribute(data->descr_A, CUSPARSE_SPMAT_DIAG_TYPE, &nonunit, sizeof(nonunit)));
        if((A.prop.uplo == 'L') == (A.order == 'R')) CHECK(cusparseSpMatSetAttribute(data->descr_A, CUSPARSE_SPMAT_FILL_MODE, &lower, sizeof(lower)));
        if((A.prop.uplo == 'U') == (A.order == 'R')) CHECK(cusparseSpMatSetAttribute(data->descr_A, CUSPARSE_SPMAT_FILL_MODE, &upper, sizeof(upper)));
    }
}



template<typename T, typename I>
void w_cusparse_trsm_dcsx_ddny_ddny<T,I>::internal_set_matrix_X()
{
    if(called_set_X) {
        CHECK(cusparseCreateDnMat(&data->descr_X, X.nrows, X.ncols, X.ld, X.vals, cusparse_data_type<T>(), cusparse_order(X.order)));
    }
    else {
        CHECK(cusparseDnMatSetValues(data->descr_X, X.vals));
    }
}



template<typename T, typename I>
void w_cusparse_trsm_dcsx_ddny_ddny<T,I>::internal_set_matrix_B()
{
    if(called_set_B) {
        CHECK(cusparseCreateDnMat(&data->descr_B, B.nrows, B.ncols, B.ld, B.vals, cusparse_data_type<T>(), cusparse_order(B.order)));
    }
    else {
        CHECK(cusparseDnMatSetValues(data->descr_B, B.vals));
    }
}



template<typename T, typename I>
void w_cusparse_trsm_dcsx_ddny_ddny<T,I>::internal_setup()
{
    data->op_A = ((A.order == 'R') ? CUSPARSE_OPERATION_NON_TRANSPOSE : CUSPARSE_OPERATION_TRANSPOSE);
    data->op_B = CUSPARSE_OPERATION_NON_TRANSPOSE;
    data->spsm_alg = CUSPARSE_SPSM_ALG_DEFAULT;

    data->handle_cusparse = dynamic_cast<handle_cusparse_new*>(handle.get())->internal.h;

    CHECK(cusparseSpSM_createDescr(&descr_spsm));

    T alpha = T{1};

    wss_internal = 0;
    CHECK(cusparseSpSM_bufferSize(data->handle_cusparse, data->op_A, data->op_B, &alpha, data->descr_A, data->descr_B, data->descr_C, cusparse_data_type<T>(), data->spsm_alg, data->descr_spsm, &wss_persistent));
    wss_tmp_preprocess = 0;
    wss_tmp_perform = 0;
}



template<typename T, typename I>
void w_cusparse_trsm_dcsx_ddny_ddny<T,I>::internal_preprocess(void * /*ws_tmp*/)
{
    T alpha = T{1};
    CHECK(cusparseSpSM_analysis(data->handle_cusparse, data->op_A, data->op_B, &alpha, data->descr_A, data->descr_B, data->descr_C, cusparse_data_type<T>(), data->spsm_alg, data->descr_spsm, ws_persistent));
}



template<typename T, typename I>
void w_cusparse_trsm_dcsx_ddny_ddny<T,I>::internal_update()
{
#if CUDART_VERSION >= 12040
    CHECK(cusparseSpSM_updateMatrix(data->handle_cusparse, data->descr_spsm, A.vals, CUSPARSE_SPSM_UPDATE_GENERAL));
#else
    T alpha = T{1};
    CHECK(cusparseSpSM_analysis(data->handle_cusparse, data->op_A, data->op_B, &alpha, data->descr_A, data->descr_B, data->descr_C, cusparse_data_type<T>(), data->spsm_alg, data->descr_spsm, ws_persistent));
#endif
}



template<typename T, typename I>
void w_cusparse_trsm_dcsx_ddny_ddny<T,I>::internal_perform(void * /*ws_tmp*/)
{
    T alpha = T{1};
    CHECK(cusparseSpSM_solve(data->handle_cusparse, data->op_A, data->op_B, &alpha, data->descr_A, data->descr_B, data->descr_C, cusparse_data_type<T>(), data->spsm_alg, data->descr_spsm));
}



#define INSTANTIATE_T_I(T,I) \
template class w_cusparse_trsm_dcsx_ddny_ddny<T,I>;

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
#endif
