
#include "gpu/operations/gemm_dcsx_ddny_ddnz.h"

#include "wrappers/cuda/operations/w_cusparse_gemm_dcsx_ddny_ddnz.h"
// #include "wrappers/rocm/operations/w_rocsparse_gemm_dcsx_ddny_ddnz.h"
// #include "wrappers/oneapi/operations/w_oneapisparse_gemm_dcsx_ddny_ddnz.h"



namespace espreso {
namespace gpu {
namespace operations {



template<typename T, typename I>
std::unique_ptr<gemm_dcsx_ddny_ddnz<T,I>> gemm_dcsx_ddny_ddnz<T,I>::make()
{
    // feel free to make this runtime ifs based on ecf or env
    #ifdef ESPRESO_USE_WRAPPER_GPU_CUDA
        return std::make_unique<w_cusparse_gemm_dcsx_ddny_ddnz<T,I>>();
    #endif
    // #ifdef ESPRESO_USE_WRAPPER_GPU_ROCM
    //     return std::make_unique<w_rocsparse_gemm_dcsx_ddny_ddnz<T,I>>();
    // #endif
    // #ifdef ESPRESO_USE_WRAPPER_GPU_ONEAPI
    //     return std::make_unique<w_oneapisparse_gemm_dcsx_ddny_ddnz<T,I>>();
    // #endif
    eslog::error("wrapper for gemm_dcsx_ddny_ddnz not available\n");
}



template<typename T, typename I>
void gemm_dcsx_ddny_ddnz<T,I>::set_handle(gpu::mgm::queue q_, gpu::spblas::handle spblas_handle_)
{
    if(called_set_handles) eslog::error("handles are already set\n");

    q = q_;
    spblas_handle = spblas_handle_;

    called_set_handles = true;
}



template<typename T, typename I>
void gemm_dcsx_ddny_ddnz<T,I>::set_matrix_A(MatrixCsxView_new<T,I> * A_)
{
    if(!called_set_handles) eslog::error("handles are not set\n");
    if(A != nullptr) eslog::error("matrix A is already set\n");
    if(A_ == nullptr) eslog::error("A cannot be nullptr\n");

    A = A_;
}



template<typename T, typename I>
void gemm_dcsx_ddny_ddnz<T,I>::set_matrix_B(MatrixDenseView_new<T> B_)
{
    if(!called_set_handles) eslog::error("handles are not set\n");
    if(B != nullptr) eslog::error("matrix B is already set\n");
    if(B_ == nullptr) eslog::error("B cannot be nullptr\n");

    B = B_;
}



template<typename T, typename I>
void gemm_dcsx_ddny_ddnz<T,I>::set_matrix_C(MatrixDenseView_new<T> C_)
{
    if(!called_set_handles) eslog::error("handles are not set\n");
    if(C != nullptr) eslog::error("matrix C is already set\n");
    if(C_ == nullptr) eslog::error("C cannot be nullptr\n");

    C = C_;
}



template<typename T, typename I>
void gemm_dcsx_ddny_ddnz<T,I>::set_coefficients(T alpha_, T beta_)
{
    alpha = alpha_;
    beta = beta_;
}



template<typename T, typename I>
void gemm_dcsx_ddny_ddnz<T,I>::setup()
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



template<typename T, typename I>
size_t gemm_dcsx_ddny_ddnz<T,I>::get_wss_internal()
{
    if(!called_setup) eslog::error("setup was not called\n");

    return wss_internal;
}



template<typename T, typename I>
size_t gemm_dcsx_ddny_ddnz<T,I>::get_wss_persistent()
{
    if(!called_setup) eslog::error("setup was not called\n");

    return wss_persistent;
}



template<typename T, typename I>
size_t gemm_dcsx_ddny_ddnz<T,I>::get_wss_tmp_preprocess()
{
    if(!called_setup) eslog::error("setup was not called\n");

    return wss_tmp_preprocess;
}



template<typename T, typename I>
size_t gemm_dcsx_ddny_ddnz<T,I>::get_wss_tmp_perform()
{
    if(!called_setup) eslog::error("setup was not called\n");

    return wss_tmp_perform;
}



template<typename T, typename I>
void gemm_dcsx_ddny_ddnz<T,I>::set_ws_persistent(void * ws_persistent_)
{
    if(ws_persistent_ == nullptr && wss_persistent > 0) eslog::error("persistent workspace is null\n");
    if(ws_persistent != nullptr) eslog::error("cannot re-set persistent workspace\n");

    ws_persistent = ws_persistent_;
}



template<typename T, typename I>
void gemm_dcsx_ddny_ddnz<T,I>::preprocess_submit(void * ws_tmp)
{
    if(!called_setup) eslog::error("setup has not been called\n");
    if(called_preprocess) eslog::error("preprocess has already been called\n");
    if(ws_tmp == nullptr && wss_tmp_preprocess > 0) eslog::error("temporary workspace is null\n");

    this->internal_preprocess(ws_tmp);

    called_preprocess = true;
}



template<typename T, typename I>
void gemm_dcsx_ddny_ddnz<T,I>::perform_submit(void * ws_tmp)
{
    if(!called_preprocess) eslog::error("preprocess has not been called\n");
    if(ws_tmp == nullptr && wss_tmp_perform > 0) eslog::error("temporary workspace is null\n");

    this->internal_perform(ws_tmp);
}



template<typename T, typename I>
void gemm_dcsx_ddny_ddnz<T,I>::submit_all(gpu::mgm::queue q, gpu::spblas::handle handle_spblas, MatrixCsxView_new<T,I> * A, MatrixDenseView_new<T> * B, MatrixDenseView_new<T> * C, T alpha, T beta, Allocator_new * ator_gpu)
{
    gemm_dcsx_ddny_ddnz<T,I> instance;
    instance.set_handles(q, handle_spblas);
    instance.set_matrix_A(A);
    instance.set_matrix_B(B);
    instance.set_matrix_C(C);
    instance.set_coefficients(alpha, beta);
    instance.setup();
    size_t wss_persistent = instance.get_wss_persistent();
    size_t wss_tmp_preprocess = instance.get_wss_tmp_preprocess();
    size_t wss_tmp_perform = instance.get_wss_tmp_perform();
    size_t wss_tmp = std::max(wss_tmp_preprocess, wss_tmp_perform);
    void * ws_persistent = ator_gpu->alloc(wss_persistent);
    void * ws_tmp = ator_gpu->alloc(wss_tmp);
    instance.set_ws_persistent(ws_persistent);
    instance.preprocess_submit(ws_tmp);
    instance.perform_submit(ws_tmp);
    gpu::mgm::submit_host_function(q, [ator_gpu,ws_persistent,ws_tmp](){
        ator_gpu->free(ws_persistent);
        ator_gpu->free(ws_tmp);
    });
}



#define INSTANTIATE_T_I(T,I) \
template class gemm_dcsx_ddny_ddnz<T,I>;

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
