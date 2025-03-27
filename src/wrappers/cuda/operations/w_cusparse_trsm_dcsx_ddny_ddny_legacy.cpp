
#ifdef HAVE_CUDA
#ifdef USE_CUSPARSE_LEGACY

#include "wrappers/cuda/operations/w_cusparse_trsm_dcsx_ddny_ddny.h"

#if defined(__GNUC__) && !defined(__clang__)
#define MY_COMPILER_GCC
#elif defined(__clang__)
#define MY_COMPILER_CLANG
#endif

// The legacy cusparse API is deprecated. I know, no need to remind me.
#ifdef MY_COMPILER_GCC
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#endif
#ifdef MY_COMPILER_CLANG
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wdeprecated-declarations"
#endif

#include "wrappers/cuda/common_cusparse.h"
#include "wrappers/cuda/common_cuda_mgm.h"



namespace espreso {
namespace gpu {
namespace operations {



struct w_cusparse_trsm_dcsx_ddny_ddny_data
{
    cusparseMatDescr_t descr_A;
    cusparseOperation_t op_A;
    cusparseOperation_t op_X;
    csrsm2Info_t info;
    cusparseSolvePolicy_t policy;
};



template<typename T, typename I>
w_cusparse_trsm_dcsx_ddny_ddny<T,I>::w_cusparse_trsm_dcsx_ddny_ddny()
{
    data = std::make_unique<w_cusparse_trsm_dcsx_ddny_ddny_data>();
}



template<typename T, typename I>
w_cusparse_trsm_dcsx_ddny_ddny<T,I>::~w_cusparse_trsm_dcsx_ddny_ddny()
{
    if(this->called_setup) {
        CHECK(cusparseDestroyMatDescr(data->descr_A));
        CHECK(cusparseDestroyCsrsm2Info(data->info));
    }

    data.reset();
}



template<typename T, typename I>
char w_cusparse_trsm_dcsx_ddny_ddny<T,I>::internal_get_native_place()
{
    return 'I';
}



template<typename T, typename I>
void w_cusparse_trsm_dcsx_ddny_ddny<T,I>::internal_setup()
{
    CHECK(cusparseCreateCsrsm2Info(&data->info));
    data->policy = CUSPARSE_SOLVE_POLICY_USE_LEVEL;
    data->op_A = ((A->order == 'R') ? CUSPARSE_OPERATION_NON_TRANSPOSE : CUSPARSE_OPERATION_TRANSPOSE);
    data->op_X = ((X->order == 'C') ? CUSPARSE_OPERATION_NON_TRANSPOSE : CUSPARSE_OPERATION_TRANSPOSE);

    auto diag = ((A->prop.diag == 'U') ? CUSPARSE_DIAG_TYPE_UNIT : CUSPARSE_DIAG_TYPE_NON_UNIT);
    auto uplo = (((A->prop.uplo == 'L') == (A->order == 'R')) ? CUSPARSE_FILL_MODE_LOWER : CUSPARSE_FILL_MODE_UPPER);

    CHECK(cusparseCreateMatDescr(&data->descr_A));
    CHECK(cusparseSetMatDiagType(data->descr_A, diag));
    CHECK(cusparseSetMatFillMode(data->descr_A, uplo));
    CHECK(cusparseSetMatIndexBase(data->descr_A, CUSPARSE_INDEX_BASE_ZERO));
    CHECK(cusparseSetMatType(data->descr_A, CUSPARSE_MATRIX_TYPE_GENERAL));

    size_t buffersize;
    using U = cpp_to_cuda_type_t<T>;
    U alpha = U{1};
    if constexpr(std::is_same_v<T,float>)                CHECK(cusparseScsrsm2_bufferSizeExt(handle_spblas->h, 1, data->op_A, data->op_X, X->nrows, X->ncols, A->nnz, &alpha, data->descr_A, (U*)A->vals, A->ptrs, A->idxs, (U*)X->vals, X->ld, data->info, data->policy, &buffersize));
    if constexpr(std::is_same_v<T,double>)               CHECK(cusparseDcsrsm2_bufferSizeExt(handle_spblas->h, 1, data->op_A, data->op_X, X->nrows, X->ncols, A->nnz, &alpha, data->descr_A, (U*)A->vals, A->ptrs, A->idxs, (U*)X->vals, X->ld, data->info, data->policy, &buffersize));
    if constexpr(std::is_same_v<T,std::complex<float>>)  CHECK(cusparseCcsrsm2_bufferSizeExt(handle_spblas->h, 1, data->op_A, data->op_X, X->nrows, X->ncols, A->nnz, &alpha, data->descr_A, (U*)A->vals, A->ptrs, A->idxs, (U*)X->vals, X->ld, data->info, data->policy, &buffersize));
    if constexpr(std::is_same_v<T,std::complex<double>>) CHECK(cusparseZcsrsm2_bufferSizeExt(handle_spblas->h, 1, data->op_A, data->op_X, X->nrows, X->ncols, A->nnz, &alpha, data->descr_A, (U*)A->vals, A->ptrs, A->idxs, (U*)X->vals, X->ld, data->info, data->policy, &buffersize));

    wss_internal = ((A->order == 'R') ? 0 : 8 * A->nnz);
    wss_persistent = 0;
    wss_tmp_preprocess = buffersize;
    wss_tmp_perform = buffersize;
}



template<typename T, typename I>
void w_cusparse_trsm_dcsx_ddny_ddny<T,I>::internal_preprocess(void * ws_tmp)
{
    using U = cpp_to_cuda_type_t<T>;
    U alpha = U{1};
    if constexpr(std::is_same_v<T,float>)                CHECK(cusparseScsrsm2_analysis(handle_spblas->h, 1, data->op_A, data->op_X, X->nrows, X->ncols, A->nnz, &alpha, data->descr_A, (U*)A->vals, A->ptrs, A->idxs, (U*)X->vals, X->ld, data->info, data->policy, ws_tmp));
    if constexpr(std::is_same_v<T,double>)               CHECK(cusparseDcsrsm2_analysis(handle_spblas->h, 1, data->op_A, data->op_X, X->nrows, X->ncols, A->nnz, &alpha, data->descr_A, (U*)A->vals, A->ptrs, A->idxs, (U*)X->vals, X->ld, data->info, data->policy, ws_tmp));
    if constexpr(std::is_same_v<T,std::complex<float>>)  CHECK(cusparseCcsrsm2_analysis(handle_spblas->h, 1, data->op_A, data->op_X, X->nrows, X->ncols, A->nnz, &alpha, data->descr_A, (U*)A->vals, A->ptrs, A->idxs, (U*)X->vals, X->ld, data->info, data->policy, ws_tmp));
    if constexpr(std::is_same_v<T,std::complex<double>>) CHECK(cusparseZcsrsm2_analysis(handle_spblas->h, 1, data->op_A, data->op_X, X->nrows, X->ncols, A->nnz, &alpha, data->descr_A, (U*)A->vals, A->ptrs, A->idxs, (U*)X->vals, X->ld, data->info, data->policy, ws_tmp));
}



template<typename T, typename I>
void w_cusparse_trsm_dcsx_ddny_ddny<T,I>::internal_perform(void * ws_tmp)
{
    if(place == 'O') {
        gpu::mgm::copy_submit(q, *B, *X);
    }

    using U = cpp_to_cuda_type_t<T>;
    U alpha = U{1};
    if constexpr(std::is_same_v<T,float>)                CHECK(cusparseScsrsm2_solve(handle_spblas->h, 1, data->op_A, data->op_X, X->nrows, X->ncols, A->nnz, &alpha, data->descr_A, (U*)A->vals, A->ptrs, A->idxs, (U*)X->vals, X->ld, data->info, data->policy, ws_tmp));
    if constexpr(std::is_same_v<T,double>)               CHECK(cusparseDcsrsm2_solve(handle_spblas->h, 1, data->op_A, data->op_X, X->nrows, X->ncols, A->nnz, &alpha, data->descr_A, (U*)A->vals, A->ptrs, A->idxs, (U*)X->vals, X->ld, data->info, data->policy, ws_tmp));
    if constexpr(std::is_same_v<T,std::complex<float>>)  CHECK(cusparseCcsrsm2_solve(handle_spblas->h, 1, data->op_A, data->op_X, X->nrows, X->ncols, A->nnz, &alpha, data->descr_A, (U*)A->vals, A->ptrs, A->idxs, (U*)X->vals, X->ld, data->info, data->policy, ws_tmp));
    if constexpr(std::is_same_v<T,std::complex<double>>) CHECK(cusparseZcsrsm2_solve(handle_spblas->h, 1, data->op_A, data->op_X, X->nrows, X->ncols, A->nnz, &alpha, data->descr_A, (U*)A->vals, A->ptrs, A->idxs, (U*)X->vals, X->ld, data->info, data->policy, ws_tmp));
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

#ifdef MY_COMPILER_GCC
#pragma GCC diagnostic pop
#endif
#ifdef MY_COMPILER_CLANG
#pragma clang diagnostic pop
#endif

#endif
#endif
