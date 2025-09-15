
#ifdef HAVE_CUDA

#include "wrappers/cuda/operations/w_cublas_hemm_ddnx_ddny_ddnz.h"

#include "wrappers/cuda/common_cublas.h"
#include "wrappers/cuda/common_cuda_mgm.h"



namespace espreso {
namespace gpu {
namespace operations {



template<typename T>
w_cublas_hemm_ddnx_ddny_ddnz<T>::w_cublas_hemm_ddnx_ddny_ddnz() = default;



template<typename T>
w_cublas_hemm_ddnx_ddny_ddnz<T>::~w_cublas_hemm_ddnx_ddny_ddnz() = default;



template<typename T>
void w_cublas_hemm_ddnx_ddny_ddnz<T>::internal_setup()
{
    // https://docs.nvidia.com/cuda/cublas/index.html#cublassetworkspace
    // no manual workspace needed if I use just a single stream with this handle, which I do
    wss_tmp_perform = 0;

    if(B->order != C->order) eslog::error("order of matrices B and C must match\n");
    if(utils::is_complex<T>() && A->order != C->order) eslog::error("for complex, order of all matrices must match\n");
}



template<typename T>
void w_cublas_hemm_ddnx_ddny_ddnz<T>::internal_perform(void * /*ws_tmp*/)
{
    cublasSideMode_t side = ((C->order == 'C') ? CUBLAS_SIDE_LEFT : CUBLAS_SIDE_RIGHT);
    cublasFillMode_t uplo = ((((C->order == 'C') == (A->order == 'C')) == (A->prop.uplo == 'L')) ? CUBLAS_FILL_MODE_LOWER : CUBLAS_FILL_MODE_UPPER);
    size_t m = C->nrows;
    size_t n = C->ncols;
    if(side == CUBLAS_SIDE_RIGHT) std::swap(m, n);

    using U = cpp_to_cuda_type_t<T>;
    if constexpr(std::is_same_v<T,float>)                CHECK(cublasSsymm(handle_dnblas->h, side, uplo, m, n, (U*)&alpha, (U*)A->vals, A->ld, (U*)B->vals, B->ld, (U*)&beta, (U*)C->vals, C->ld));
    if constexpr(std::is_same_v<T,double>)               CHECK(cublasDsymm(handle_dnblas->h, side, uplo, m, n, (U*)&alpha, (U*)A->vals, A->ld, (U*)B->vals, B->ld, (U*)&beta, (U*)C->vals, C->ld));
    if constexpr(std::is_same_v<T,std::complex<float>>)  CHECK(cublasChemm(handle_dnblas->h, side, uplo, m, n, (U*)&alpha, (U*)A->vals, A->ld, (U*)B->vals, B->ld, (U*)&beta, (U*)C->vals, C->ld));
    if constexpr(std::is_same_v<T,std::complex<double>>) CHECK(cublasZhemm(handle_dnblas->h, side, uplo, m, n, (U*)&alpha, (U*)A->vals, A->ld, (U*)B->vals, B->ld, (U*)&beta, (U*)C->vals, C->ld));
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
