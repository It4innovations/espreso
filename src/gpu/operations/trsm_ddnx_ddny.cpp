
#include "gpu/operations/trsm_ddnx_ddny.h"

#include "basis/utilities/stacktimer.h"
#include "wrappers/cuda/operations/w_cublas_trsm_ddnx_ddny.h"



namespace espreso {
namespace gpu {
namespace operations {



template<typename T>
std::unique_ptr<trsm_ddnx_ddny<T>> trsm_ddnx_ddny<T>::make()
{
    #ifdef ESPRESO_USE_WRAPPER_GPU_CUDA
        return std::make_unique<w_cublas_trsm_ddnx_ddny<T>>();
    #endif
    eslog::error("wrapper for trsm_ddnx_ddny not available\n");
}



template<typename T>
void trsm_ddnx_ddny<T>::set_handles(gpu::mgm::queue q_, gpu::dnblas::handle handle_dnblas_)
{
    if(called_set_handles) eslog::error("handles are already set\n");

    q = q_;
    handle_dnblas = handle_dnblas_;

    called_set_handles = true;
}



template<typename T>
void trsm_ddnx_ddny<T>::set_matrix_A(MatrixDenseView_new<T> * A_)
{
    if(!called_set_handles) eslog::error("handles are not set\n");
    if(A != nullptr) eslog::error("matrix A is already set\n");
    if(A_ == nullptr) eslog::error("A cannot be nullptr\n");

    A = A_;
}



template<typename T>
void trsm_ddnx_ddny<T>::set_matrix_X(MatrixDenseView_new<T> * X_)
{
    if(!called_set_handles) eslog::error("handles are not set\n");
    if(X != nullptr) eslog::error("matrix X is already set\n");
    if(X_ == nullptr) eslog::error("X cannot be nullptr\n");

    X = X_;
}



template<typename T>
void trsm_ddnx_ddny<T>::setup()
{
    stacktimer::push("trsm_ddnx_ddny::setup");

    if(!called_set_handles) eslog::error("handles are not set\n");
    if(A == nullptr) eslog::error("matrix A is not set\n");
    if(X == nullptr) eslog::error("matrix X is not set\n");
    if(called_setup) eslog::error("setup was already called\n");
    if(A->nrows != A->ncols) eslog::error("matrix A is not square\n");
    if(A->nrows != X->nrows) eslog::error("incompatible matrices\n");
    if(A->prop.uplo != 'L' && A->prop.uplo != 'U') eslog::error("invalid matrix A uplo\n");

    this->internal_setup();

    stacktimer::pop();

    called_setup = true;
}



template<typename T>
size_t trsm_ddnx_ddny<T>::get_wss_tmp_perform()
{
    if(!called_setup) eslog::error("setup was not called\n");

    return wss_tmp_perform;
}



template<typename T>
void trsm_ddnx_ddny<T>::perform_submit(void * ws_tmp)
{
    stacktimer::push("trsm_ddnx_ddny::perform_submit");

    if(!called_setup) eslog::error("setup was not called\n");
    if(ws_tmp == nullptr && wss_tmp_perform > 0) eslog::error("temporary workspace is null\n");

    this->internal_perform(ws_tmp);

    stacktimer::pop();
}



#define INSTANTIATE_T(T) \
template class trsm_ddnx_ddny<T>;

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
