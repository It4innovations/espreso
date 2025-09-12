
#include "gpu/operations/convert_dcsx_ddny.h"

#include "basis/utilities/stacktimer.h"
#include "wrappers/cuda/operations/w_cusparse_convert_dcsx_ddny.h"
#include "wrappers/rocm/operations/w_rocsparse_convert_dcsx_ddny.h"



namespace espreso {
namespace gpu {
namespace operations {



template<typename T, typename I>
std::unique_ptr<convert_dcsx_ddny<T,I>> convert_dcsx_ddny<T,I>::make()
{
    // feel free to make this runtime ifs based on ecf or env
    #ifdef ESPRESO_USE_WRAPPER_GPU_CUDA
        return std::make_unique<w_cusparse_convert_dcsx_ddny<T,I>>();
    #endif
    #ifdef ESPRESO_USE_WRAPPER_GPU_ROCM
        return std::make_unique<w_rocsparse_convert_dcsx_ddny<T,I>>();
    #endif
    eslog::error("wrapper for convert_dcsx_ddny not available\n");
}



template<typename T, typename I>
void convert_dcsx_ddny<T,I>::set_handles(gpu::mgm::queue q_, gpu::spblas::handle handle_spblas_)
{
    if(called_set_handles) eslog::error("handles are already set\n");

    q = q_;
    handle_spblas = handle_spblas_;

    called_set_handles = true;
}



template<typename T, typename I>
void convert_dcsx_ddny<T,I>::set_matrix_src(MatrixCsxView_new<T,I> * M_src_)
{
    if(!called_set_handles) eslog::error("handles are not set\n");
    if(M_src != nullptr) eslog::error("source matrix is already set\n");
    if(M_src_ == nullptr) eslog::error("source matrix cannot be nullptr\n");

    M_src = M_src_;
}



template<typename T, typename I>
void convert_dcsx_ddny<T,I>::set_matrix_dst(MatrixDenseView_new<T> * M_dst_)
{
    if(!called_set_handles) eslog::error("handles are not set\n");
    if(M_dst != nullptr) eslog::error("destination matrix is already set\n");
    if(M_dst_ == nullptr) eslog::error("destination matrix cannot be nullptr\n");

    M_dst = M_dst_;
}



template<typename T, typename I>
void convert_dcsx_ddny<T,I>::setup()
{
    stacktimer::push("convert_dcsx_ddny::setup");

    if(!called_set_handles) eslog::error("handles are not set\n");
    if(M_src == nullptr) eslog::error("source matrix is not set\n");
    if(M_dst == nullptr) eslog::error("destination matrix is not set\n");
    if(!M_src->ator->is_data_accessible_gpu()) eslog::error("source matrix must be gpu-accessible\n");
    if(!M_dst->ator->is_data_accessible_gpu()) eslog::error("destination matrix must be gpu-accessible\n");
    if(called_setup) eslog::error("setup was already called\n");
    if(M_src->nrows != M_dst->nrows || M_src->ncols != M_dst->ncols) eslog::error("matrix sizes dont match\n");

    this->internal_setup();

    // stacktimer::info("wss_internal       %zu", wss_internal);
    // stacktimer::info("wss_persistent     %zu", wss_persistent);
    // stacktimer::info("wss_tmp_preprocess %zu", wss_tmp_preprocess);
    // stacktimer::info("wss_tmp_perform    %zu", wss_tmp_perform);

    stacktimer::pop();

    called_setup = true;
}



template<typename T, typename I>
size_t convert_dcsx_ddny<T,I>::get_wss_internal()
{
    if(!called_setup) eslog::error("setup was not called\n");

    return wss_internal;
}



template<typename T, typename I>
size_t convert_dcsx_ddny<T,I>::get_wss_persistent()
{
    if(!called_setup) eslog::error("setup was not called\n");

    return wss_persistent;
}



template<typename T, typename I>
size_t convert_dcsx_ddny<T,I>::get_wss_tmp_preprocess()
{
    if(!called_setup) eslog::error("setup was not called\n");

    return wss_tmp_preprocess;
}



template<typename T, typename I>
size_t convert_dcsx_ddny<T,I>::get_wss_tmp_perform()
{
    if(!called_setup) eslog::error("setup was not called\n");

    return wss_tmp_perform;
}



template<typename T, typename I>
void convert_dcsx_ddny<T,I>::set_ws_persistent(void * ws_persistent_)
{
    if(ws_persistent_ == nullptr && wss_persistent > 0) eslog::error("persistent workspace is null\n");
    if(ws_persistent != nullptr) eslog::error("cannot re-set persistent workspace\n");

    ws_persistent = ws_persistent_;
}



template<typename T, typename I>
void convert_dcsx_ddny<T,I>::preprocess_submit(void * ws_tmp)
{
    stacktimer::push("convert_dcsx_ddny::preprocess_submit");

    if(!called_setup) eslog::error("setup has not been called\n");
    if(called_preprocess) eslog::error("preprocess has already been called\n");
    if(ws_tmp == nullptr && wss_tmp_preprocess > 0) eslog::error("temporary workspace is null\n");

    this->internal_preprocess(ws_tmp);

    stacktimer::pop();

    called_preprocess = true;
}



template<typename T, typename I>
void convert_dcsx_ddny<T,I>::perform_submit(void * ws_tmp)
{
    stacktimer::push("convert_dcsx_ddny::perform_submit");

    if(!called_preprocess) eslog::error("preprocess has not been called\n");
    if(ws_tmp == nullptr && wss_tmp_perform > 0) eslog::error("temporary workspace is null\n");

    this->internal_perform(ws_tmp);

    stacktimer::pop();
}



#define INSTANTIATE_T_I(T,I) \
template class convert_dcsx_ddny<T,I>;

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
