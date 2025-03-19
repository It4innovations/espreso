
#include "gpu/operations/copy_ddnx_ddnx.h"

#include "wrappers/cuda/operations/w_cuda_copy_ddnx_ddnx.h"



namespace espreso {
namespace gpu {
namespace operations {



template<typename T, typename I>
std::unique_ptr<copy_ddnx_ddnx<T,I>> copy_ddnx_ddnx<T,I>::make()
{
    #ifdef ESPRESO_USE_WRAPPER_GPU_CUDA
        return std::make_unique<w_cuda_copy_ddnx_ddnx<T>>();
    #endif
    eslog::error("wrapper for copy_ddnx_ddnx not available\n");
}



template<typename T, typename I>
void copy_ddnx_ddnx<T,I>::set_handles(gpu::mgm::queue q_)
{
    if(called_set_handles) eslog::error("handles are already set\n");

    q = q_;

    called_set_handles = true;
}



template<typename T, typename I>
void copy_ddnx_ddnx<T,I>::set_matrix_src(MatrixDenseView_new<T,I> * M_src_)
{
    if(M_src != nullptr) eslog::error("matrix M_src is already set\n");
    if(M_src_ == nullptr) eslog::error("M_src cannot be nullptr\n");

    M_src = M_src_;
}



template<typename T, typename I>
void copy_ddnx_ddnx<T,I>::set_matrix_dst(MatrixDenseView_new<T,I> * M_dst_)
{
    if(M_dst != nullptr) eslog::error("matrix M_dst is already set\n");
    if(M_dst_ == nullptr) eslog::error("M_dst cannot be nullptr\n");

    M_dst = M_dst_;
}



template<typename T, typename I>
void copy_ddnx_ddnx<T,I>::set_uplo(char uplo_);
{
    uplo = uplo_;
}



template<typename T, typename I>
void copy_ddnx_ddnx<T,I>::setup()
{
    if(!called_set_handles) eslog::error("handles are not set\n");
    if(M_src == nullptr) eslog::error("source matrix is not sen\n");
    if(M_dst == nullptr) eslog::error("destination matrix is not set\n");
    if(M_src->nrows != M_dst->nrows || M_src->ncols != M_dst->ncols) eslog::error("matrix sizes dont match\n");
    if(M_src->order != M_dst->order) eslog::error("matrix orders must match\n");
    if(M_src->prop.uplo != M_dst->prop.uplo) eslog::error("matrix uplo must match\n");

    this->internal_setup();

    called_setup = true;
}



template<typename T, typename I>
size_t copy_ddnx_ddnx<T,I>::get_wss_tmp_perform()
{
    if(!called_setup) eslog::error("setup was not called\n");

    return wss_tmp_perform;
}



template<typename T, typename I>
void copy_ddnx_ddnx<T,I>::perform_submit(void * ws_tmp)
{
    if(!called_setup) eslog::error("setup was not called\n");
    if(ws_tmp == nullptr && wss_tmp_perform > 0) eslog::error("temporary workspace is null\n");

    this->internal_perform(ws_tmp);
}




#define INSTANTIATE_T(T) \
template class copy_ddnx_ddnx<T>;

    #define INSTANTIATE \
    /* INSTANTIATE_T(float) */ \
    INSTANTIATE_T(double) \
    /* INSTANTIATE_T(std::complex<float>) */ \
    /* INSTANTIATE_T(std::complex<double>) */

        INSTANTIATE

    #undef INSTANTIATE
#undef INSTANTIATE_T


    
}
}
}
