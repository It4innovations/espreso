
#include "gpu/operations/hemm_ddnx_ddny_ddnz.h"

#include "basis/utilities/stacktimer.h"
#include "wrappers/cuda/operations/w_cublas_hemm_ddnx_ddny_ddnz.h"



namespace espreso {
namespace gpu {
namespace operations {



template<typename T>
std::unique_ptr<hemm_ddnx_ddny_ddnz<T>> hemm_ddnx_ddny_ddnz<T>::make()
{
    // feel free to make this runtime ifs based on ecf or env
    #ifdef ESPRESO_USE_WRAPPER_GPU_CUDA
        return std::make_unique<w_cublas_hemm_ddnx_ddny_ddnz<T>>();
    #endif
    eslog::error("wrapper for hemm_ddnx_ddny_ddnz not available\n");
}



template<typename T>
void hemm_ddnx_ddny_ddnz<T>::set_handles(gpu::mgm::queue q_, gpu::dnblas::handle handle_dnblas_)
{
    if(called_set_handles) eslog::error("handles are already set\n");

    q = q_;
    handle_dnblas = handle_dnblas_;

    called_set_handles = true;
}



template<typename T>
void hemm_ddnx_ddny_ddnz<T>::set_matrix_A(MatrixDenseView_new<T> * A_)
{
    if(A != nullptr) eslog::error("matrix A is already set\n");
    if(A_ == nullptr) eslog::error("A cannot be nullptr\n");

    A = A_;
}



template<typename T>
void hemm_ddnx_ddny_ddnz<T>::set_matrix_B(MatrixDenseView_new<T> * B_)
{
    if(B != nullptr) eslog::error("matrix B is already set\n");
    if(B_ == nullptr) eslog::error("B cannot be nullptr\n");

    B = B_;
}



template<typename T>
void hemm_ddnx_ddny_ddnz<T>::set_matrix_C(MatrixDenseView_new<T> * C_)
{
    if(C != nullptr) eslog::error("matrix C is already set\n");
    if(C_ == nullptr) eslog::error("C cannot be nullptr\n");

    C = C_;
}



template<typename T>
void hemm_ddnx_ddny_ddnz<T>::set_coefficients(T alpha_, T beta_)
{
    alpha = alpha_;
    beta = beta_;
}



template<typename T>
void hemm_ddnx_ddny_ddnz<T>::setup()
{
    stacktimer::push("hemm_ddnx_ddny_ddnz::setup");

    if(!called_set_handles) eslog::error("handles are not set\n");
    if(A == nullptr) eslog::error("matrix A is not set\n");
    if(B == nullptr) eslog::error("matrix B is not set\n");
    if(C == nullptr) eslog::error("matrix C is not set\n");
    if(!A->ator->is_data_accessible_gpu()) eslog::error("matrix A must be gpu-accessible\n");
    if(!B->ator->is_data_accessible_gpu()) eslog::error("matrix B must be gpu-accessible\n");
    if(!C->ator->is_data_accessible_gpu()) eslog::error("matrix C must be gpu-accessible\n");
    if(called_setup) eslog::error("setup was already called\n");
    if(!is_hermitian<T>(A->prop.symm)) eslog::error("matrix A must be hermitian\n");
    if(A->nrows != C->nrows || B->ncols != C->ncols || A->ncols != B->nrows) eslog::error("incompatible matrices\n");

    this->internal_setup();

    // stacktimer::info("wss_tmp_perform    %zu", wss_tmp_perform);

    stacktimer::pop();

    called_setup = true;
}



template<typename T>
size_t hemm_ddnx_ddny_ddnz<T>::get_wss_tmp_perform()
{
    if(!called_setup) eslog::error("setup was not called\n");

    return wss_tmp_perform;
}



template<typename T>
void hemm_ddnx_ddny_ddnz<T>::perform_submit(void * ws_tmp)
{
    stacktimer::push("hemm_ddnx_ddny_ddnz::perform_submit");

    if(!called_setup) eslog::error("setup has not been called\n");
    if(ws_tmp == nullptr && wss_tmp_perform > 0) eslog::error("temporary workspace is null\n");

    this->internal_perform(ws_tmp);

    stacktimer::pop();
}



#define INSTANTIATE_T(T) \
template class hemm_ddnx_ddny_ddnz<T>;

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
