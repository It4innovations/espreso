
#include "gpu/operations/convert_ddnx_ddny.h"

#include "basis/utilities/stacktimer.h"
#include "wrappers/cuda/operations/w_cuda_convert_ddnx_ddny.h"



namespace espreso {
namespace gpu {
namespace operations {



template<typename T>
std::unique_ptr<convert_ddnx_ddny<T>> convert_ddnx_ddny<T>::make()
{
    #ifdef ESPRESO_USE_WRAPPER_GPU_CUDA
        return std::make_unique<w_cuda_convert_ddnx_ddny<T>>();
    #endif
    eslog::error("wrapper for convert_ddnx_ddny not available\n");
}



template<typename T>
void convert_ddnx_ddny<T>::set_handles(gpu::mgm::queue q_, gpu::dnblas::handle handle_dnblas_)
{
    if(called_set_handles) eslog::error("handles are already set\n");

    q = q_;
    handle_dnblas = handle_dnblas_;

    called_set_handles = true;
}



template<typename T>
void convert_ddnx_ddny<T>::set_matrix_src(MatrixDenseView_new<T> * M_src_)
{
    if(M_src != nullptr) eslog::error("matrix M_src is already set\n");
    if(M_src_ == nullptr) eslog::error("M_src cannot be nullptr\n");

    M_src = M_src_;
}



template<typename T>
void convert_ddnx_ddny<T>::set_matrix_dst(MatrixDenseView_new<T> * M_dst_)
{
    if(M_dst != nullptr) eslog::error("matrix M_dst is already set\n");
    if(M_dst_ == nullptr) eslog::error("M_dst cannot be nullptr\n");

    M_dst = M_dst_;
}



template<typename T>
void convert_ddnx_ddny<T>::setup()
{
    stacktimer::push("convert_ddnx_ddny::setup");

    if(!called_set_handles) eslog::error("handles are not set\n");
    if(M_src == nullptr) eslog::error("source matrix is not sen\n");
    if(M_dst == nullptr) eslog::error("destination matrix is not set\n");
    if(M_src->nrows != M_dst->nrows || M_src->ncols != M_dst->ncols) eslog::error("matrix sizes dont match\n");

    this->internal_setup();

    stacktimer::pop();

    called_setup = true;
}



template<typename T>
size_t convert_ddnx_ddny<T>::get_wss_tmp_perform()
{
    if(!called_setup) eslog::error("setup was not called\n");

    return wss_tmp_perform;
}



template<typename T>
void convert_ddnx_ddny<T>::perform_submit(void * ws_tmp)
{
    stacktimer::push("convert_ddnx_ddny::perform_submit");

    if(!called_setup) eslog::error("setup was not called\n");
    if(ws_tmp == nullptr && wss_tmp_perform > 0) eslog::error("temporary workspace is null\n");

    this->internal_perform(ws_tmp);

    stacktimer::pop();
}




#define INSTANTIATE_T(T) \
template class convert_ddnx_ddny<T>;

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
