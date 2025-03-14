
#include "gpu/operations/convert_dcsx_ddny.h"

#include "wrappers/cuda/operations/w_cusparse_convert_dcsx_ddny.h"



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
    eslog::error("wrapper for convert_dcsx_ddny not available\n");
}



template<typename T, typename I>
void convert_dcsx_ddny<T,I>::set_handles(gpu::mgm::queue q_, gpu::spblas::handle handle_spblas_)
{
    if(called_set_handles) eslog::error("handles are already set\n");

    q = q_;
    spblas_handle = spblas_handle_;

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
    if(!called_set_handles) eslog::error("handles are not set\n");
    if(M_src == nullptr) eslog::error("source matrix is not set\n");
    if(M_dst == nullptr) eslog::error("destination matrix is not set\n");
    if(called_setup) eslog::error("setup was already called\n");
    if(M_src->nrows != M_dst->nrows || M_src->ncols != M_dst->ncols) eslog::error("matrix sizes dont match\n");

    this->internal_setup();

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
    if(!called_setup) eslog::error("setup has not been called\n");
    if(called_preprocess) eslog::error("preprocess has already been called\n");
    if(ws_tmp == nullptr && wss_tmp_preprocess > 0) eslog::error("temporary workspace is null\n");

    this->internal_preprocess(ws_tmp);

    called_preprocess = true;
}



template<typename T, typename I>
void convert_dcsx_ddny<T,I>::perform_submit(void * ws_tmp)
{
    if(!called_preprocess) eslog::error("preprocess has not been called\n");
    if(ws_tmp == nullptr && wss_tmp_perform > 0) eslog::error("temporary workspace is null\n");

    this->internal_perform(ws_tmp);
}



template<typename T, typename I>
void convert_dcsx_ddny<T,I>::submit_all(gpu::mgm::queue q, gpu::spblas::handle handle_spblas, MatrixCsxView_new<T,I> * M_src, MatrixDenseView_new<T> * M_dst, Allocator_new * ator_gpu)
{
    auto instance = convert_dcsx_ddny<T,I>::make();
    instance->set_handles(q, handle_spblas);
    instance->set_matrix_src(M_src);
    instance->set_matrix_dst(M_dst);
    instance->setup();
    size_t wss_persistent = instance->get_wss_persistent();
    size_t wss_tmp_preprocess = instance->get_wss_tmp_preprocess();
    size_t wss_tmp_perform = instance->get_wss_tmp_perform();
    size_t wss_tmp = std::max(wss_tmp_preprocess, wss_tmp_perform);
    void * ws_persistent = nullptr;
    void * ws_tmp = nullptr;
    if(wss_persistent > 0) ws_persistent = ator_gpu->alloc(wss_persistent);
    if(wss_tmp > 0) ws_tmp = ator_gpu->alloc(wss_tmp);
    instance->set_ws_persistent(ws_persistent);
    insatnce->preprocess_submit(ws_tmp);
    instance->perform_submit(ws_tmp);
    if(ws_persistent != nullptr || ws_tmp != nullptr) {
        gpu::mgm::submit_host_function(q, [ator_gpu,ws_persistent,ws_tmp](){
            ator_gpu->free(ws_persistent);
            ator_gpu->free(ws_tmp);
        });
    }
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



#endif /* SRC_GPU_OPERATIONS_CONVERT_CSX_DNY_H */
