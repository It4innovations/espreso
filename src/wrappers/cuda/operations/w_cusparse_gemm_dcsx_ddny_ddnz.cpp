
#ifdef HAVE_CUDA

#include "wrappers/cuda/operations/w_cusparse_gemm_dcsx_ddny_ddnz.h"

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



template<typename T, typename I>
w_cusparse_gemm_dcsx_ddny_ddnz<T,I>::w_cusparse_gemm_dcsx_ddny_ddnz() = default;



template<typename T, typename I>
w_cusparse_gemm_dcsx_ddny_ddnz<T,I>::~w_cusparse_gemm_dcsx_ddny_ddnz()
{
    if(this->called_setup) {
        CHECK(cusparseDestroySpMat(data->descr_A));
        CHECK(cusparseDestroyDnMat(data->descr_B));
        CHECK(cusparseDestroyDnMat(data->descr_C));
    }
}



template<typename T, typename I>
void w_cusparse_gemm_dcsx_ddny_ddnz<T,I>::internal_setup()
{
    data = std::make_unique<w_cusparse_gemm_dcsx_ddny_ddnz_data<T,I>>();

    data->op_A = CUSPARSE_OPERATION_NON_TRANSPOSE;
    data->op_B = CUSPARSE_OPERATION_NON_TRANSPOSE;
    data->spmm_alg = CUSPARSE_SPMM_ALG_DEFAULT;

    data->handle_cusparse = handle_spblas->h;

    T * dummyptr_T = (T*)(sizeof(T));
    I * dummyptr_I = (I*)(sizeof(I));
    if(A->order == 'R') CHECK(cusparseCreateCsr(&data->descr_A, A->nrows, A->ncols, A->nnz, dummyptr_I, dummyptr_I, dummyptr_T, cusparse_index_type<I>(), cusparse_index_type<I>(), CUSPARSE_INDEX_BASE_ZERO, cusparse_data_type<T>()));
    if(A->order == 'C') CHECK(cusparseCreateCsc(&data->descr_A, A->nrows, A->ncols, A->nnz, dummyptr_I, dummyptr_I, dummyptr_T, cusparse_index_type<I>(), cusparse_index_type<I>(), CUSPARSE_INDEX_BASE_ZERO, cusparse_data_type<T>()));
    CHECK(cusparseCreateDnMat(&data->descr_B, B->nrows, B->ncols, B->ld, dummyptr_T, cusparse_data_type<T>(), cusparse_order(B->order)));
    CHECK(cusparseCreateDnMat(&data->descr_C, C->nrows, C->ncols, C->ld, dummyptr_T, cusparse_data_type<T>(), cusparse_order(C->order)));

    size_t buffersize;
    CHECK(cusparseSpMM_bufferSize(data->handle_cusparse, data->op_A, data->op_B, &alpha, data->descr_A, data->descr_B, &beta, data->descr_C, cusparse_data_type<T>(), data->spmm_alg, &buffersize));

    wss_internal = 0; // TODO: check
    wss_persistent = utils::round_up(buffersize, gpu::mgm::get_natural_pitch_align());
    wss_tmp_preprocess = 0;
    wss_tmp_perform = 0;
}



template<typename T, typename I>
void w_cusparse_gemm_dcsx_ddny_ddnz<T,I>::internal_preprocess(void * /*ws_tmp*/)
{
    if(A->order == 'R') CHECK(cusparseCsrSetPointers(data->descr_A, A->ptrs, A->idxs, A->vals));
    if(A->order == 'C') CHECK(cusparseCscSetPointers(data->descr_A, A->ptrs, A->idxs, A->vals));

    CHECK(cusparseSpMM_preprocess(data->handle_cusparse, data->op_A, data->op_B, &alpha, data->descr_A, data->descr_B, &beta, data->descr_C, cusparse_data_type<T>(), data->spmm_alg, ws_persistent));
}



template<typename T, typename I>
void w_cusparse_gemm_dcsx_ddny_ddnz<T,I>::internal_perform(void * /*ws_tmp*/)
{
    if(A->order == 'R') CHECK(cusparseCsrSetPointers(data->descr_A, A->ptrs, A->idxs, A->vals));
    if(A->order == 'C') CHECK(cusparseCscSetPointers(data->descr_A, A->ptrs, A->idxs, A->vals));
    CHECK(cusparseDnMatSetValues(data->descr_B, B->vals));
    CHECK(cusparseDnMatSetValues(data->descr_C, C->vals));

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
