
#include "gpu/operations/herk_ddnx_ddny.h"

#include "wrappers/cuda/operations/w_cublas_herk_ddnx_ddny.h"



namespace espreso {
namespace gpu {
namespace operations {



template<typename T>
std::unique_ptr<herk_ddnx_ddny<T>> herk_ddnx_ddny<T>::make()
{
    #ifdef ESPRESO_USE_WRAPPER_GPU_CUDA
        return std::make_unique<w_cublas_herk_ddnx_ddny<T>>();
    #endif
    eslog::error("wrapper for herk_ddnx_ddny not available\n");
}



template<typename T>
void herk_ddnx_ddny<T>::set_handles(gpu::mgm::queue q_, gpu::dnblas::handle dnblas_handle_)
{
    if(called_set_handles) eslog::error("handles are already set\n");

    q = q_;
    dnblas_handle = dnblas_handle_;

    called_set_handles = true;
}



template<typename T>
void herk_ddnx_ddny<T>::set_matrix_A(MatrixDenseView_new<T,I> A_)
{
    if(!called_set_handles) eslog::error("handles are not set\n");
    if(called_set_A && !MatrixDenseView_new<T>::are_interchangable(A, A_)) eslog::error("invalid replacement for matrix A\n");

    A = A_;

    internal_set_matrix_A();

    called_set_A = true;
}



template<typename T>
void herk_ddnx_ddny<T>::set_matrix_C(MatrixDenseView_new<T> C_)
{
    if(!called_set_handles) eslog::error("handles are not set\n");
    if(called_set_C && !MatrixDenseView_new<T>::are_interchangable(C, C_)) eslog::error("invalid replacement for matrix C\n");

    C = C_;

    internal_set_matrix_C();

    called_set_C = true;
}



template<typename T>
void herk_ddnx_ddny<T>::set_coefficients(Treal alpha_, Treal beta_)
{
    alpha = alpha_;
    beta = beta_;
}



template<typename T>
void herk_ddnx_ddny<T>::set_mode(herk_mode mode_)
{
    mode = mode_;

    called_set_mode = true;
}



template<typename T>
void herk_ddnx_ddny<T>::setup()
{
    if(!called_set_handles) eslog::error("handles are not set\n");
    if(!called_set_A) eslog::error("matrix A is not set\n");
    if(!called_set_C) eslog::error("matrix C is not set\n");
    if(!called_set_mode) eslog::error("mode is not set\n");
    if(called_setup) eslog::error("setup was already called\n");
    if(C.nrows != C.ncols) eslog::error("matrix C is not square\n");
    if(mode == herk_mode::AhA && A.ncols != C.ncols) eslog::error("incompatible matrices\n");
    if(mode == herk_mode::AAh && A.nrows != C.nrows) eslog::error("incompatible matrices\n");
    if(C.prop.uplo != 'L' && C.prop.uplo != 'U') eslog::error("invalid matrix A uplo\n");

    this->internal_setup();

    called_setup = true;
}



template<typename T>
size_t herk_ddnx_ddny<T>::get_wss_tmp_perform()
{
    if(!called_setup) eslog::error("setup was not called\n");

    return wss_tmp_perform;
}



template<typename T>
void herk_ddnx_ddny<T>::perform_submit(void * ws_tmp)
{
    if(!called_setup) eslog::error("setup was not called\n");
    if(ws_tmp == nullptr && wss_tmp_perform > 0) eslog::error("temporary workspace is null\n");

    this->internal_perform(ws_tmp);
}



template<typename T>
void herk_ddnx_ddny<T>::submit_all(gpu::mgm::queue q, gpu::dnblas::handle handle_dnblas, MatrixDenseView_new<T> A, MatrixDenseView_new<T> C, Treal alpha, Treal beta, herk_mode mode, Allocator_new * ator_gpu)
{
    trsm_ddnx_ddny<T> instance;
    instance.set_handles(q, handle_dnblas);
    instance.set_matrix_A(A);
    instance.set_matrix_C(C);
    instance.set_coefficients(alpha, beta);
    instance.set_mode(mode);
    instance.setup();
    size_t wss_tmp = instance.get_wss_tmp_perform();
    void * ws_tmp = nullptr;
    if(wss_tmp > 0) {
        ator_gpu->alloc(wss_tmp);
    }
    instance.perform_submit(ws_tmp);
    if(wss_tmp > 0) {
        gpu::mgm::submit_host_function(q, [ator_gpu,ws_tmp](){
            ator_gpu->free(ws_tmp);
        });
    }
}




#define INSTANTIATE_T(T) \
template class herk_ddnx_ddny<T>;

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
