
#ifdef HAVE_CUDA

#include "wrappers/cuda/operations/w_cublas_herk_ddnx_ddny.h"

#include "wrappers/cuda/common_cublas.h"
#include "wrappers/cuda/common_cuda_mgm.h"



namespace espreso {
namespace gpu {
namespace operations {



struct w_cublas_herk_ddnx_ddny_data
{
    cublasFillMode_t uplo_C;
    cublasOperation_t trans_mode;
    size_t n;
    size_t k;
};



template<typename T>
w_cublas_herk_ddnx_ddny<T>::w_cublas_herk_ddnx_ddny()
{
    data = std::make_unique<w_cublas_herk_ddnx_ddny_data>();
}



template<typename T>
w_cublas_herk_ddnx_ddny<T>::~w_cublas_herk_ddnx_ddny()
{
    data.reset();
}



template<typename T>
void w_cublas_herk_ddnx_ddny<T>::internal_setup()
{
    // https://docs.nvidia.com/cuda/cublas/index.html#cublassetworkspace
    // no manual workspace needed if I use just a single stream with this handle, which I do
    wss_tmp_perform = 0;

    data->uplo_C = (((C->prop.uplo == 'L') == (C->order == 'C')) ? CUBLAS_FILL_MODE_LOWER : CUBLAS_FILL_MODE_UPPER);
    data->trans_mode = (((mode == math::blas::herk_mode::AAh) == (A->order == 'C')) ? CUBLAS_OP_N : CUBLAS_OP_T);
    data->n = A->nrows;
    data->k = A->ncols;
    if(mode == math::blas::herk_mode::AhA) {
        std::swap(data->n, data->k);
    }
}



template<typename T>
void w_cublas_herk_ddnx_ddny<T>::internal_perform(void * /*ws_tmp*/)
{
    if(utils::is_complex<T>()) eslog::error("complex types not yet supported in herk\n");
    using U = cpp_to_cuda_type_t<T>;
    if constexpr(std::is_same_v<T,float>)                CHECK(cublasSsyrk(handle_dnblas->h, data->uplo_C, data->trans_mode, data->n, data->k, &alpha, (U*)A->vals, A->ld, &beta, (U*)C->vals, C->ld));
    if constexpr(std::is_same_v<T,double>)               CHECK(cublasDsyrk(handle_dnblas->h, data->uplo_C, data->trans_mode, data->n, data->k, &alpha, (U*)A->vals, A->ld, &beta, (U*)C->vals, C->ld));
    if constexpr(std::is_same_v<T,std::complex<float>>)  CHECK(cublasCherk(handle_dnblas->h, data->uplo_C, data->trans_mode, data->n, data->k, &alpha, (U*)A->vals, A->ld, &beta, (U*)C->vals, C->ld));
    if constexpr(std::is_same_v<T,std::complex<double>>) CHECK(cublasZherk(handle_dnblas->h, data->uplo_C, data->trans_mode, data->n, data->k, &alpha, (U*)A->vals, A->ld, &beta, (U*)C->vals, C->ld));
}



#define INSTANTIATE_T(T) \
template class w_cublas_herk_ddnx_ddny<T>;

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
