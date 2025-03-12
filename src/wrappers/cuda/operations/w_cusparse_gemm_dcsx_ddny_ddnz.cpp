
#ifdef HAVE_CUDA

#include "wrappers/cuda/operations/w_cusparse_gemm_dcsx_ddny_ddnz.h"

#include <cusparse.h>

#include "wrappers/cuda/common_cusparse.h"



namespace espreso {
namespace gpu {
namespace operations {



template<typename T, typename I>
struct w_cusparse_gemm_dcsx_ddny_ddnz_data
{
    cusparseHandle_t handle_cusparse;
    cusparseSpMatDescr_t descr_A;
    cusparseDnMatDescr_t descr_B;
    cusparseDnMatDescr_t descr_C;
    cusparseOperation_t op_A;
    cusparseOperation_t op_B;
    cusparseSpMMAlg_t spmm_alg;
};



template<typename T, typaname I>
w_cusparse_gemm_dcsx_ddny_ddnz<T,I>::~w_cusparse_gemm_dcsx_ddny_ddnz()
{
    data = std::make_unique<w_cusparse_gemm_dcsx_ddny_ddnz_data<T,I>>();
}



template<typename T, typaname I>
w_cusparse_gemm_dcsx_ddny_ddnz<T,I>::~w_cusparse_gemm_dcsx_ddny_ddnz()
{
    if(called_set_A) {
        CHECK(cusparseDestroySpMat(data->descr_A));
    }
    if(called_set_B) {
        CHECK(cusparseDestroyDnMat(data->descr_B));
    }
    if(called_set_C) {
        CHECK(cusparseDestroyDnMat(data->descr_C));
    }

    data->reset();
}



template<typename T, typaname I>
void w_cusparse_gemm_dcsx_ddny_ddnz<T,I>::internal_set_matrix_A()
{
    if(called_set_A) {
        if(A.order == 'R') CHECK(cusparseCsrSetPointers(data->descr_A, A.ptrs, A.idxs, A.vals));
        if(A.order == 'C') CHECK(cusparseCscSetPointers(data->descr_A, A.ptrs, A.idxs, A.vals));
    }
    else {
        if(A.order == 'R') CHECK(cusparseCreateCsr(data->descr_A, A.nrows, A.ncols, A.nnz, A.ptrs, A.idxs, A.vals, cusparse_index_type<I>(), cusparse_index_type<I>(), CUSPARSE_INDEX_BASE_ZERO, cusparse_data_type<T>()));
        if(A.order == 'C') CHECK(cusparseCreateCsc(data->descr_A, A.nrows, A.ncols, A.nnz, A.ptrs, A.idxs, A.vals, cusparse_index_type<I>(), cusparse_index_type<I>(), CUSPARSE_INDEX_BASE_ZERO, cusparse_data_type<T>()));
    }
}



template<typename T, typaname I>
void w_cusparse_gemm_dcsx_ddny_ddnz<T,I>::internal_set_matrix_B()
{
    if(called_set_B) {
        CHECK(cusparseCreateDnMat(&data->descr_B, B.nrows, B.ncols, B.ld, B.vals, cusparse_data_type<T>(), cusparse_order(B.order)));
    }
    else {
        CHECK(cusparseDnMatSetValues(data->descr_B, B.vals));
    }
}



template<typename T, typaname I>
void w_cusparse_gemm_dcsx_ddny_ddnz<T,I>::internal_set_matrix_C()
{
    if(called_set_C) {
        CHECK(cusparseCreateDnMat(&data->descr_C, C.nrows, C.ncols, C.ld, C.vals, cusparse_data_type<T>(), cusparse_order(C.order)));
    }
    else {
        CHECK(cusparseDnMatSetValues(data->descr_C, C.vals));
    }
}



template<typename T, typaname I>
void w_cusparse_gemm_dcsx_ddny_ddnz<T,I>::internal_setup()
{
    data->op_A = CUSPARSE_OPERATION_NON_TRANSPOSE;
    data->op_B = CUSPARSE_OPERATION_NON_TRANSPOSE;
    data->spmm_alg = CUSPARSE_SPMM_ALG_DEFAULT;

    data->handle_cusparse = handle_spblas->h;

    wss_internal = 0; // TODO check
    CHECK(cusparseSpMM_bufferSize(data->handle_cusparse, data->op_A, data->op_B, &alpha, data->descr_A, data->descr_B, &beta, data->descr_C, cusparse_data_type<T>(), data->spmm_alg, &wss_persistent));
    wss_tmp_preprocess = 0;
    wss_tmp_perform = 0;
}



template<typename T, typaname I>
void w_cusparse_gemm_dcsx_ddny_ddnz<T,I>::internal_preprocess(void * /*ws_tmp*/)
{
    CHECK(cusparseSpMM_preprocess(data->handle_cusparse, data->op_A, data->op_B, &alpha, data->descr_A, data->descr_B, &beta, data->descr_C, cusparse_data_type<T>(), data->spmm_alg, ws_persistent));
}



template<typename T, typaname I>
void w_cusparse_gemm_dcsx_ddny_ddnz<T,I>::internal_perform(void * /*ws_tmp*/)
{
    CHECK(cusparseSpMM(data->handle_cusparse, data->op_A, data->op_B, &alpha, data->descr_A, data->descr_B, &beta, data->descr_C, cusparse_data_type<T>(), data->spmm_alg, ws_persistent));
}



#define INSTANTIATE_T_I(T,I) \
template class w_cusparse_gemm_dcsx_ddny_ddnz<T,I>;

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
