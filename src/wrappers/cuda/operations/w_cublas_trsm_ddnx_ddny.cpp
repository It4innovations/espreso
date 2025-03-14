
#ifdef HAVE_CUDA

#include "wrappers/cuda/operations/w_cublas_trsm_ddnx_ddny.h"



namespace espreso {
namespace gpu {
namespace operations {



struct w_cublas_trsm_ddnx_ddny_data
{
    cublasSideMode_t side;
    cublasFillMode_t uplo_A;
    cublasOperation_t op_A;
    cublasDiagType_t diag_A;
    size_t m;
    size_t n;
};



template<typename T>
w_cublas_trsm_ddnx_ddny<T>::w_cublas_trsm_ddnx_ddny()
{
    data = std::make_unique<w_cublas_trsm_ddnx_ddny_data>();
}



template<typename T>
w_cublas_trsm_ddnx_ddny<T>::~w_cublas_trsm_ddnx_ddny()
{
    data.reset();
}



template<typename T>
void w_cublas_trsm_ddnx_ddny<T>::internal_setup()
{
    // https://docs.nvidia.com/cuda/cublas/index.html#cublassetworkspace
    // no manual workspace needed if I use just a single stream with this handle, which I do
    wss_tmp_perform = 0;

    data->side = ((X->order == 'C') ? CUBLAS_SIDE_LEFT : CUBLAS_SIDE_RIGHT);
    data->op_A = ((X->order == A->order) ? CUBLAS_OP_N : CUBLAS_OP_T);
    data->uplo_A = (((A->prop.uplo == 'L') == (A->order == 'C')) ? CUBLAS_FILL_MODE_LOWER : CUBLAS_FILL_MODE_UPPER);
    data->m = X->nrows;
    data->n = X->ncols;
    if(X->order == 'R') {
        std::swap(data->m, data->n);
    }
}



template<typename T>
void w_cublas_trsm_ddnx_ddny<T>::internal_perform(void * ws_tmp)
{
    using U = cpp_to_cublas_type_t<T>;
    U one = U{1};
    if constexpr(std::is_same_v<T,float>)                CHECK(cublasStrsm(handle_dnblas->h, data->side, data->uplo_A, data->op_A, data->diag_A, data->m, data->n, &one, (U*)A->vals, A->ld, (U*)X->vals, X->ld));
    if constexpr(std::is_same_v<T,double>)               CHECK(cublasDtrsm(handle_dnblas->h, data->side, data->uplo_A, data->op_A, data->diag_A, data->m, data->n, &one, (U*)A->vals, A->ld, (U*)X->vals, X->ld));
    if constexpr(std::is_same_v<T,std::complex<float>>)  CHECK(cublasCtrsm(handle_dnblas->h, data->side, data->uplo_A, data->op_A, data->diag_A, data->m, data->n, &one, (U*)A->vals, A->ld, (U*)X->vals, X->ld));
    if constexpr(std::is_same_v<T,std::complex<double>>) CHECK(cublasZtrsm(handle_dnblas->h, data->side, data->uplo_A, data->op_A, data->diag_A, data->m, data->n, &one, (U*)A->vals, A->ld, (U*)X->vals, X->ld));
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
