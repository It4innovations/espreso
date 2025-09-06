
#ifdef HAVE_CUDA

#include "wrappers/cuda/operations/w_cublas_hemm_ddnx_ddny_ddnz.h"

#include "wrappers/cuda/common_cublas.h"
#include "wrappers/cuda/common_cuda_mgm.h"



namespace espreso {
namespace gpu {
namespace operations {



struct w_cublas_hemm_ddnx_ddny_ddnz_data
{
    cublasSideMode_t side;
    cublasFillMode_t uplo;
    size_t m;
    size_t n;
};



template<typename T>
w_cublas_hemm_ddnx_ddny_ddnz<T>::w_cublas_hemm_ddnx_ddny_ddnz()
{
    data = std::make_unique<w_cublas_hemm_ddnx_ddny_ddnz_data>();
}



template<typename T>
w_cublas_hemm_ddnx_ddny_ddnz<T>::~w_cublas_hemm_ddnx_ddny_ddnz()
{
    data.reset();
}



template<typename T>
void w_cublas_hemm_ddnx_ddny_ddnz<T>::internal_setup()
{
    // https://docs.nvidia.com/cuda/cublas/index.html#cublassetworkspace
    // no manual workspace needed if I use just a single stream with this handle, which I do
    wss_tmp_perform = 0;

    if(B->order != C->order) eslog::error("order of all matrices must match\n");
    if(utils::is_complex<T>() && A->order != C->order) eslog::error("for complex, order of all matrices must match\n");

    data->side = ((C->order == 'C') ? CUBLAS_SIDE_LEFT : CUBLAS_SIDE_RIGHT);
    data->uplo = ((((C->order == 'C') == (A->order == 'C')) == (A->prop.uplo == 'L')) ? CUBLAS_FILL_MODE_LOWER : CUBLAS_FILL_MODE_UPPER);

    data->m = C->nrows;
    data->n = C->ncols;
    if(data->side == CUBLAS_SIDE_RIGHT) std::swap(data->m, data->n);
}



template<typename T>
void w_cublas_hemm_ddnx_ddny_ddnz<T>::internal_perform(void * /*ws_tmp*/)
{
    using U = cpp_to_cuda_type_t<T>;
    if constexpr(std::is_same_v<T,float>)                CHECK(cublasSsymm(handle_dnblas->h, data->side, data->uplo, data->m, data->n, (U*)&alpha, (U*)A->vals, A->ld, (U*)B->vals, B->ld, (U*)&beta, (U*)C->vals, C->ld));
    if constexpr(std::is_same_v<T,double>)               CHECK(cublasDsymm(handle_dnblas->h, data->side, data->uplo, data->m, data->n, (U*)&alpha, (U*)A->vals, A->ld, (U*)B->vals, B->ld, (U*)&beta, (U*)C->vals, C->ld));
    if constexpr(std::is_same_v<T,std::complex<float>>)  CHECK(cublasChemm(handle_dnblas->h, data->side, data->uplo, data->m, data->n, (U*)&alpha, (U*)A->vals, A->ld, (U*)B->vals, B->ld, (U*)&beta, (U*)C->vals, C->ld));
    if constexpr(std::is_same_v<T,std::complex<double>>) CHECK(cublasZhemm(handle_dnblas->h, data->side, data->uplo, data->m, data->n, (U*)&alpha, (U*)A->vals, A->ld, (U*)B->vals, B->ld, (U*)&beta, (U*)C->vals, C->ld));
}



#define INSTANTIATE_T(T) \
template class w_cublas_hemm_ddnx_ddny_ddnz<T>;

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
