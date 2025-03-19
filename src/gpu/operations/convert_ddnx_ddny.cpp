
#include "gpu/operations/convert_ddnx_ddny.h"

#include "wrappers/cuda/operations/w_cuda_convert_ddnx_ddny.h"



namespace espreso {
namespace gpu {
namespace operations {



template<typename T, typename I>
std::unique_ptr<convert_ddnx_ddny<T,I>> convert_ddnx_ddny<T,I>::make()
{
    #ifdef ESPRESO_USE_WRAPPER_GPU_CUDA
        return std::make_unique<w_cuda_convert_ddnx_ddny<T>>();
    #endif
    eslog::error("wrapper for convert_ddnx_ddny not available\n");
}



template<typename T, typename I>
void convert_ddnx_ddny<T,I>::set_handles(gpu::mgm::queue q_, gpu::dnblas::handle handle_dnblas_)
{
    if(called_set_handles) eslog::error("handles are already set\n");

    q = q_;
    handle_dnblas = handle_dnblas_;

    called_set_handles = true;
}



template<typename T, typename I>
void convert_ddnx_ddny<T,I>::set_matrix_src(MatrixDenseView_new<T,I> * M_src_)
{
    if(M_src != nullptr) eslog::error("matrix M_src is already set\n");
    if(M_src_ == nullptr) eslog::error("M_src cannot be nullptr\n");

    M_src = M_src_;
}



template<typename T, typename I>
void convert_ddnx_ddny<T,I>::set_matrix_dst(MatrixDenseView_new<T,I> * M_dst_)
{
    if(M_dst != nullptr) eslog::error("matrix M_dst is already set\n");
    if(M_dst_ == nullptr) eslog::error("M_dst cannot be nullptr\n");

    M_dst = M_dst_;
}



template<typename T, typename I>
void convert_ddnx_ddny<T,I>::setup()
{
    if(!called_set_handles) eslog::error("handles are not set\n");
    if(M_src == nullptr) eslog::error("source matrix is not sen\n");
    if(M_dst == nullptr) eslog::error("destination matrix is not set\n");
    if(M_src->nrows != M_dst->nrows || M_src->ncols != M_dst->ncols) eslog::error("matrix sizes dont match\n");

    this->internal_setup();

    called_setup = true;
}



template<typename T, typename I>
size_t convert_ddnx_ddny<T,I>::get_wss_tmp_perform()
{
    if(!called_setup) eslog::error("setup was not called\n");

    return wss_tmp_perform;
}



template<typename T, typename I>
void convert_ddnx_ddny<T,I>::perform_submit(void * ws_tmp)
{
    if(!called_setup) eslog::error("setup was not called\n");
    if(ws_tmp == nullptr && wss_tmp_perform > 0) eslog::error("temporary workspace is null\n");

    this->internal_perform(ws_tmp);
}




#define INSTANTIATE_T(T) \
template class handle_dnblas<T>;

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
