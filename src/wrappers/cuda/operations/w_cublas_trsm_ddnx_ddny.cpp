
#ifdef HAVE_CUDA

#include "wrappers/cuda/operations/w_cublas_trsm_ddnx_ddny.h"

#include "wrappers/cuda/common_cublas.h"
#include "wrappers/cuda/common_cuda_mgm.h"



namespace espreso {
namespace gpu {
namespace operations {



template<typename T>
w_cublas_trsm_ddnx_ddny<T>::w_cublas_trsm_ddnx_ddny() = default;



template<typename T>
w_cublas_trsm_ddnx_ddny<T>::~w_cublas_trsm_ddnx_ddny() = default;



template<typename T>
void w_cublas_trsm_ddnx_ddny<T>::internal_setup()
{
    // https://docs.nvidia.com/cuda/cublas/index.html#cublassetworkspace
    // no manual workspace needed if I use just a single stream with this handle, which I do
    wss_tmp_perform = 0;
}



template<typename T>
void w_cublas_trsm_ddnx_ddny<T>::internal_perform(void * ws_tmp)
{
    cublasSideMode_t side = ((X->order == 'C') ? CUBLAS_SIDE_LEFT : CUBLAS_SIDE_RIGHT);
    cublasOperation_t op_A = ((X->order == A->order) ? CUBLAS_OP_N : CUBLAS_OP_T);
    cublasFillMode_t uplo_A = (((A->prop.uplo == 'L') == (A->order == 'C')) ? CUBLAS_FILL_MODE_LOWER : CUBLAS_FILL_MODE_UPPER);
    cublasDiagType_t diag_A = ((A->prop.diag == 'U') ? CUBLAS_DIAG_UNIT : CUBLAS_DIAG_NON_UNIT);
    size_t m = X->nrows;
    size_t n = X->ncols;
    if(X->order == 'R') {
        std::swap(m, n);
    }

    using U = cpp_to_cuda_type_t<T>;
    U one = U{1};
    if constexpr(std::is_same_v<T,float>)                CHECK(cublasStrsm(handle_dnblas->h, side, uplo_A, op_A, diag_A, m, n, &one, (U*)A->vals, A->ld, (U*)X->vals, X->ld));
    if constexpr(std::is_same_v<T,double>)               CHECK(cublasDtrsm(handle_dnblas->h, side, uplo_A, op_A, diag_A, m, n, &one, (U*)A->vals, A->ld, (U*)X->vals, X->ld));
    if constexpr(std::is_same_v<T,std::complex<float>>)  CHECK(cublasCtrsm(handle_dnblas->h, side, uplo_A, op_A, diag_A, m, n, &one, (U*)A->vals, A->ld, (U*)X->vals, X->ld));
    if constexpr(std::is_same_v<T,std::complex<double>>) CHECK(cublasZtrsm(handle_dnblas->h, side, uplo_A, op_A, diag_A, m, n, &one, (U*)A->vals, A->ld, (U*)X->vals, X->ld));
}



#define INSTANTIATE_T(T) \
template class w_cublas_trsm_ddnx_ddny<T>;

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
