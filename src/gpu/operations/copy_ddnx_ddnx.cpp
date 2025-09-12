
#include "gpu/operations/copy_ddnx_ddnx.h"

#include "basis/utilities/stacktimer.h"
#include "wrappers/cuda/operations/w_cuda_copy_ddnx_ddnx.h"
#include "wrappers/rocm/operations/w_rocm_copy_ddnx_ddnx.h"



namespace espreso {
namespace gpu {
namespace operations {



template<typename T>
std::unique_ptr<copy_ddnx_ddnx<T>> copy_ddnx_ddnx<T>::make()
{
    #ifdef ESPRESO_USE_WRAPPER_GPU_CUDA
        return std::make_unique<w_cuda_copy_ddnx_ddnx<T>>();
    #endif
    #ifdef ESPRESO_USE_WRAPPER_GPU_ROCM
        return std::make_unique<w_rocm_copy_ddnx_ddnx<T>>();
    #endif
    eslog::error("wrapper for copy_ddnx_ddnx not available\n");
}



template<typename T>
void copy_ddnx_ddnx<T>::set_handles(gpu::mgm::queue q_)
{
    if(called_set_handles) eslog::error("handles are already set\n");

    q = q_;

    called_set_handles = true;
}



template<typename T>
void copy_ddnx_ddnx<T>::set_matrix_src(MatrixDenseView_new<T> * M_src_)
{
    if(M_src != nullptr) eslog::error("matrix M_src is already set\n");
    if(M_src_ == nullptr) eslog::error("M_src cannot be nullptr\n");

    M_src = M_src_;
}



template<typename T>
void copy_ddnx_ddnx<T>::set_matrix_dst(MatrixDenseView_new<T> * M_dst_)
{
    if(M_dst != nullptr) eslog::error("matrix M_dst is already set\n");
    if(M_dst_ == nullptr) eslog::error("M_dst cannot be nullptr\n");

    M_dst = M_dst_;
}



template<typename T>
void copy_ddnx_ddnx<T>::set_uplo(char uplo_)
{
    uplo = uplo_;
}



template<typename T>
void copy_ddnx_ddnx<T>::setup()
{
    stacktimer::push("copy_ddnx_ddnx::setup");

    if(!called_set_handles) eslog::error("handles are not set\n");
    if(M_src == nullptr) eslog::error("source matrix is not sen\n");
    if(M_dst == nullptr) eslog::error("destination matrix is not set\n");
    if(!M_src->ator->is_data_accessible_gpu()) eslog::error("source matrix must be gpu-accessible\n");
    if(!M_dst->ator->is_data_accessible_gpu()) eslog::error("destination matrix must be gpu-accessible\n");
    if(M_src->nrows != M_dst->nrows || M_src->ncols != M_dst->ncols) eslog::error("matrix sizes dont match\n");
    if(M_src->order != M_dst->order) eslog::error("matrix orders must match\n");

    this->internal_setup();

    // stacktimer::info("wss_tmp_perform    %zu", wss_tmp_perform);

    stacktimer::pop();

    called_setup = true;
}



template<typename T>
size_t copy_ddnx_ddnx<T>::get_wss_tmp_perform()
{
    if(!called_setup) eslog::error("setup was not called\n");

    return wss_tmp_perform;
}



template<typename T>
void copy_ddnx_ddnx<T>::perform_submit(void * ws_tmp)
{
    stacktimer::push("copy_ddnx_ddnx::perform_submit");

    if(!called_setup) eslog::error("setup was not called\n");
    if(ws_tmp == nullptr && wss_tmp_perform > 0) eslog::error("temporary workspace is null\n");

    this->internal_perform(ws_tmp);

    stacktimer::pop();
}




#define INSTANTIATE_T(T) \
template class copy_ddnx_ddnx<T>;

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
