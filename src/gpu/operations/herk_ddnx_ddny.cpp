
#include "gpu/operations/herk_ddnx_ddny.h"

#include "basis/utilities/stacktimer.h"
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
void herk_ddnx_ddny<T>::set_handles(gpu::mgm::queue q_, gpu::dnblas::handle handle_dnblas_)
{
    if(called_set_handles) eslog::error("handles are already set\n");

    q = q_;
    handle_dnblas = handle_dnblas_;

    called_set_handles = true;
}



template<typename T>
void herk_ddnx_ddny<T>::set_matrix_A(MatrixDenseView_new<T> * A_)
{
    if(A != nullptr) eslog::error("matrix A is already set\n");
    if(A_ == nullptr) eslog::error("A cannot be nullptr\n");

    A = A_;
}



template<typename T>
void herk_ddnx_ddny<T>::set_matrix_C(MatrixDenseView_new<T> * C_)
{
    if(C != nullptr) eslog::error("matrix C is already set\n");
    if(C_ == nullptr) eslog::error("C cannot be nullptr\n");

    C = C_;
}



template<typename T>
void herk_ddnx_ddny<T>::set_coefficients(Treal alpha_, Treal beta_)
{
    alpha = alpha_;
    beta = beta_;
}



template<typename T>
void herk_ddnx_ddny<T>::set_mode(math::blas::herk_mode mode_)
{
    mode = mode_;

    called_set_mode = true;
}



template<typename T>
void herk_ddnx_ddny<T>::setup()
{
    stacktimer::push("herk_ddnx_ddny::setup");

    if(!called_set_handles) eslog::error("handles are not set\n");
    if(A == nullptr) eslog::error("matrix A is not set\n");
    if(C == nullptr) eslog::error("matrix C is not set\n");
    if(!called_set_mode) eslog::error("mode is not set\n");
    if(called_setup) eslog::error("setup was already called\n");
    if(C->nrows != C->ncols) eslog::error("matrix C is not square\n");
    if(mode == math::blas::herk_mode::AhA && A->ncols != C->ncols) eslog::error("incompatible matrices\n");
    if(mode == math::blas::herk_mode::AAh && A->nrows != C->nrows) eslog::error("incompatible matrices\n");
    if(C->prop.uplo != 'L' && C->prop.uplo != 'U') eslog::error("invalid matrix C uplo\n");

    this->internal_setup();

    // stacktimer::info("wss_tmp_perform    %zu", wss_tmp_perform);

    stacktimer::pop();

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
    stacktimer::push("herk_ddnx_ddny::perform_submit");

    if(!called_setup) eslog::error("setup was not called\n");
    if(ws_tmp == nullptr && wss_tmp_perform > 0) eslog::error("temporary workspace is null\n");

    this->internal_perform(ws_tmp);

    stacktimer::pop();
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
