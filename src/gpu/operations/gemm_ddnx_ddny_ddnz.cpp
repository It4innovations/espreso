
#include "gpu/operations/gemm_ddnx_ddny_ddnz.h"

#include "wrappers/cuda/operations/w_cublas_gemm_ddnx_ddny_ddnz.h"



namespace espreso {
namespace gpu {
namespace operations {



template<typename T>
std::unique_ptr<gemm_ddnx_ddny_ddnz<T>> gemm_ddnx_ddny_ddnz<T>::make()
{
    #ifdef ESPRESO_USE_WRAPPER_GPU_CUDA
        return std::make_unique<w_cublas_gemm_ddnx_ddny_ddnz<T>>();
    #endif
    eslog::error("wrapper for gemm_ddnx_ddny_ddnz not available\n");
}



template<typename T>
void gemm_ddnx_ddny_ddnz<T>::set_handles(gpu::mgm::queue q_, gpu::dnblas::handle handle_dnblas_)
{
    if(called_set_handles) eslog::error("handles are already set\n");

    q = q_;
    handle_dnblas = handle_dnblas_;

    called_set_handles = true;
}



template<typename T>
void gemm_ddnx_ddny_ddnz<T>::set_matrix_A(MatrixDenseView_new<T> * A_)
{
    if(!called_set_handles) eslog::error("handles are not set\n");
    if(A != nullptr) eslog::error("matrix A is already set\n");
    if(A_ == nullptr) eslog::error("A cannot be nullptr\n");

    A = A_;
}



template<typename T>
void gemm_ddnx_ddny_ddnz<T>::set_matrix_B(MatrixDenseView_new<T> * B_)
{
    if(!called_set_handles) eslog::error("handles are not set\n");
    if(B != nullptr) eslog::error("matrix B is already set\n");
    if(B_ == nullptr) eslog::error("B cannot be nullptr\n");

    B = B_;
}



template<typename T>
void gemm_ddnx_ddny_ddnz<T>::set_matrix_C(MatrixDenseView_new<T> * C_)
{
    if(!called_set_handles) eslog::error("handles are not set\n");
    if(C != nullptr) eslog::error("matrix C is already set\n");
    if(C_ == nullptr) eslog::error("C cannot be nullptr\n");

    C = C_;
}



template<typename T>
void gemm_ddnx_ddny_ddnz<T>::set_coefficients(T alpha_, T beta_)
{
    alpha = alpha_;
    beta = beta_;
}



template<typename T>
void gemm_ddnx_ddny_ddnz<T>::setup()
{
    if(!called_set_handles) eslog::error("handles are not set\n");
    if(A == nullptr) eslog::error("matrix A is not set\n");
    if(B == nullptr) eslog::error("matrix B is not set\n");
    if(C == nullptr) eslog::error("matrix C is not set\n");
    if(called_setup) eslog::error("setup was already called\n");
    if(A->nrows != C->nrows || B->ncols != C->ncols || A->ncols != B->nrows) eslog::error("incompatible matrices\n");

    this->internal_setup();

    called_setup = true;
}



template<typename T>
size_t gemm_ddnx_ddny_ddnz<T>::get_wss_tmp_perform()
{
    if(!called_setup) eslog::error("setup was not called\n");

    return wss_tmp_perform;
}



template<typename T>
void gemm_ddnx_ddny_ddnz<T>::perform_submit(void * ws_tmp)
{
    if(!called_setup) eslog::error("setup was not called\n");
    if(ws_tmp == nullptr && wss_tmp_perform > 0) eslog::error("temporary workspace is null\n");

    this->internal_perform(ws_tmp);
}




#define INSTANTIATE_T(T) \
template class gemm_ddnx_ddny_ddnz<T>;

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
