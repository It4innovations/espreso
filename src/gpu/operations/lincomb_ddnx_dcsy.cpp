
#include "gpu/operations/lincomb_ddnx_dcsy.h"

#include "basis/utilities/stacktimer.h"
#include "wrappers/cuda/operations/w_cuda_lincomb_ddnx_dcsy.h"
#include "wrappers/rocm/operations/w_rocm_lincomb_ddnx_dcsy.h"
#include "wrappers/oneapi/operations/w_oneapi_lincomb_ddnx_dcsy.h"



namespace espreso {
namespace gpu {
namespace operations {



template<typename T, typename I>
std::unique_ptr<lincomb_ddnx_dcsy<T,I>> lincomb_ddnx_dcsy<T,I>::make()
{
    // feel free to make this runtime ifs based on ecf or env
    #ifdef ESPRESO_USE_WRAPPER_GPU_CUDA
        return std::make_unique<w_cuda_lincomb_ddnx_dcsy<T,I>>();
    #endif
    #ifdef ESPRESO_USE_WRAPPER_GPU_ROCM
        return std::make_unique<w_rocm_lincomb_ddnx_dcsy<T,I>>();
    #endif
    #ifdef ESPRESO_USE_WRAPPER_GPU_ONEAPI
        return std::make_unique<w_oneapi_lincomb_ddnx_dcsy<T,I>>();
    #endif
    eslog::error("wrapper for lincomb_ddnx_dcsy not available\n");
}



template<typename T, typename I>
void lincomb_ddnx_dcsy<T,I>::set_handles(gpu::mgm::queue q_)
{
    if(called_set_handles) eslog::error("handles are already set\n");

    q = q_;

    called_set_handles = true;
}



template<typename T, typename I>
void lincomb_ddnx_dcsy<T,I>::set_matrix_X(MatrixDenseView_new<T> * X_)
{
    if(X != nullptr) eslog::error("matrix X is already set\n");

    X = X_;
}



template<typename T, typename I>
void lincomb_ddnx_dcsy<T,I>::set_matrix_A(MatrixDenseView_new<T> * A_)
{
    if(A != nullptr) eslog::error("matrix A is already set\n");

    A = A_;
}



template<typename T, typename I>
void lincomb_ddnx_dcsy<T,I>::set_matrix_B(MatrixCsxView_new<T,I> * B_)
{
    if(B != nullptr) eslog::error("matrix B is already set\n");

    B = B_;
}



template<typename T, typename I>
void lincomb_ddnx_dcsy<T,I>::set_coefficients(Treal alpha_, Treal beta_)
{
    alpha = alpha_;
    beta = beta_;
}



template<typename T, typename I>
void lincomb_ddnx_dcsy<T,I>::setup()
{
    stacktimer::push("lincomb_ddnx_dcsy::setup");

    if(!called_set_handles) eslog::error("handles are not set\n");
    if(called_setup) eslog::error("setup was already called\n");
    if(X == nullptr) eslog::error("matrix X is not set\n");
    if(A == nullptr) eslog::error("matrix A is not set\n");
    if(B == nullptr) eslog::error("matrix B is not set\n");
    if(!X->ator->is_data_accessible_gpu()) eslog::error("matrix X must be gpu-accessible\n");
    if(!A->ator->is_data_accessible_gpu()) eslog::error("matrix A must be gpu-accessible\n");
    if(!B->ator->is_data_accessible_gpu()) eslog::error("matrix B must be gpu-accessible\n");
    if(A->nrows != X->nrows || A->ncols != X->ncols) eslog::error("non-matching size X vs A\n");
    if(B->nrows != X->nrows || B->ncols != X->ncols) eslog::error("non-matching size X vs B\n");
    if(A->prop.uplo != X->prop.uplo || B->prop.uplo != X->prop.uplo) eslog::error("non-matching uplo\n");

    this->internal_setup();

    // stacktimer::info("wss_tmp_perform    %zu", wss_tmp_perform);

    stacktimer::pop();

    called_setup = true;
}



template<typename T, typename I>
size_t lincomb_ddnx_dcsy<T,I>::get_wss_tmp_perform()
{
    return wss_tmp_perform;
}



template<typename T, typename I>
void lincomb_ddnx_dcsy<T,I>::perform_submit(void * ws_tmp)
{
    stacktimer::push("lincomb_ddnx_dcsy::perform_submit");

    if(!called_setup) eslog::error("setup has not been called\n");
    if(ws_tmp == nullptr && wss_tmp_perform > 0) eslog::error("temporary workspace is null\n");

    this->internal_perform(ws_tmp);

    stacktimer::pop();
}



#define INSTANTIATE_T_I(T,I) \
template class lincomb_ddnx_dcsy<T,I>;

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
