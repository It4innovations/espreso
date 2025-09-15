
#ifdef HAVE_CUDA

#include "wrappers/cuda/operations/w_cublas_gemm_ddnx_ddny_ddnz.h"

#include "wrappers/cuda/common_cublas.h"
#include "wrappers/cuda/common_cuda_mgm.h"



namespace espreso {
namespace gpu {
namespace operations {



template<typename T>
w_cublas_gemm_ddnx_ddny_ddnz<T>::w_cublas_gemm_ddnx_ddny_ddnz() = default;



template<typename T>
w_cublas_gemm_ddnx_ddny_ddnz<T>::~w_cublas_gemm_ddnx_ddny_ddnz() = default;



template<typename T>
void w_cublas_gemm_ddnx_ddny_ddnz<T>::internal_setup()
{
    // https://docs.nvidia.com/cuda/cublas/index.html#cublassetworkspace
    // no manual workspace needed if I use just a single stream with this handle, which I do
    wss_tmp_perform = 0;
}



template<typename T>
void w_cublas_gemm_ddnx_ddny_ddnz<T>::internal_perform(void * /*ws_tmp*/)
{
    bool swap_a_b = (C->order == 'R');
    cublasOperation_t op_A = (((A->order == 'C') == (C->order == 'C')) ? CUBLAS_OP_N : CUBLAS_OP_T);
    cublasOperation_t op_B = (((B->order == 'C') == (C->order == 'C')) ? CUBLAS_OP_N : CUBLAS_OP_T);
    size_t m = C->nrows;
    size_t n = C->ncols;
    size_t k = A->ncols;
    if(swap_a_b) {
        std::swap(m, n);
        std::swap(op_A, op_B);
    }

    using U = cpp_to_cuda_type_t<T>;
    U * A_vals = (U*)A->vals;
    U * B_vals = (U*)B->vals;
    U * C_vals = (U*)C->vals;
    size_t A_ld = A->ld;
    size_t B_ld = B->ld;
    size_t C_ld = C->ld;
    if(swap_a_b) {
        std::swap(A_vals, B_vals);
        std::swap(A_ld, B_ld);
    }

    if constexpr(std::is_same_v<T,float>)                CHECK(cublasSgemm(handle_dnblas->h, op_A, op_B, m, n, k, (U*)&alpha, A_vals, A_ld, B_vals, B_ld, (U*)&beta, C_vals, C_ld));
    if constexpr(std::is_same_v<T,double>)               CHECK(cublasDgemm(handle_dnblas->h, op_A, op_B, m, n, k, (U*)&alpha, A_vals, A_ld, B_vals, B_ld, (U*)&beta, C_vals, C_ld));
    if constexpr(std::is_same_v<T,std::complex<float>>)  CHECK(cublasCgemm(handle_dnblas->h, op_A, op_B, m, n, k, (U*)&alpha, A_vals, A_ld, B_vals, B_ld, (U*)&beta, C_vals, C_ld));
    if constexpr(std::is_same_v<T,std::complex<double>>) CHECK(cublasZgemm(handle_dnblas->h, op_A, op_B, m, n, k, (U*)&alpha, A_vals, A_ld, B_vals, B_ld, (U*)&beta, C_vals, C_ld));
}



#define INSTANTIATE_T(T) \
template class w_cublas_gemm_ddnx_ddny_ddnz<T>;

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
