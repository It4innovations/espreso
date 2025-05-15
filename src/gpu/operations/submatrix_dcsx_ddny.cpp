
#include "gpu/operations/submatrix_dcsx_ddny.h"

#include "basis/utilities/stacktimer.h"
#include "wrappers/cuda/operations/w_cuda_submatrix_dcsx_ddny.h"



namespace espreso {
namespace gpu {
namespace operations {



template<typename T, typename I>
std::unique_ptr<submatrix_dcsx_ddny<T,I>> submatrix_dcsx_ddny<T,I>::make()
{
    // feel free to make this runtime ifs based on ecf or env
    #ifdef ESPRESO_USE_WRAPPER_GPU_CUDA
        return std::make_unique<w_cuda_submatrix_dcsx_ddny<T,I>>();
    #endif
    eslog::error("wrapper for submatrix_dcsx_dcsx not available\n");
}



template<typename T, typename I>
void submatrix_dcsx_ddny<T,I>::set_handles(gpu::mgm::queue q_)
{
    if(called_set_handles) eslog::error("handles are already set\n");

    q = q_;

    called_set_handles = true;
}



template<typename T, typename I>
void submatrix_dcsx_ddny<T,I>::set_bounds(size_t row_start_, size_t row_end_, size_t col_start_, size_t col_end_)
{
    if(called_set_bounds) eslog::error("bounds are already set\n");

    row_start = row_start_;
    row_end = row_end_;
    col_start = col_start_;
    col_end = col_end_;

    called_set_bounds = true;
}



template<typename T, typename I>
void submatrix_dcsx_ddny<T,I>::set_matrix_src(MatrixCsxView_new<T,I> * M_src_)
{
    if(M_src != nullptr) eslog::error("source matrix is already set\n");
    if(M_src_ == nullptr) eslog::error("source matrix cannot be nullptr\n");

    M_src = M_src_;
}



template<typename T, typename I>
void submatrix_dcsx_ddny<T,I>::set_matrix_dst(MatrixDenseView_new<T> * M_dst_)
{
    if(M_dst != nullptr) eslog::error("destination matrix is already set\n");
    if(M_dst_ == nullptr) eslog::error("destination matrix cannot be nullptr\n");

    M_dst = M_dst_;
}



template<typename T, typename I>
void submatrix_dcsx_ddny<T,I>::setup()
{
    stacktimer::push("submatrix_dcsx_ddny::setup");

    if(!called_set_handles) eslog::error("handles are not set\n");
    if(!called_set_bounds) eslog::error("bounds are not set\n");
    if(M_src == nullptr) eslog::error("source matrix is not set\n");
    if(M_dst == nullptr) eslog::error("destination matrix is not set\n");
    if(!M_src->ator->is_data_accessible_gpu()) eslog::error("source matrix must be gpu-accessible\n");
    if(!M_dst->ator->is_data_accessible_gpu()) eslog::error("destination matrix must be gpu-accessible\n");
    if(called_setup) eslog::error("setup was already called\n");
    if(row_start > row_end || row_end > M_src->nrows || col_start > col_end || col_end > M_src->ncols) eslog::error("wrong bounds\n");
    if((row_end - row_start) != M_dst->nrows || (col_end - col_start) != M_dst->ncols) eslog::error("wrong output matrix size\n");

    primary_start = ((M_src->order == 'R') ? row_start : col_start);
    primary_end = ((M_src->order == 'R') ? row_end : col_end);
    secdary_start = ((M_src->order == 'R') ? col_start : row_start);
    secdary_end = ((M_src->order == 'R') ? col_end : row_end);

    this->internal_setup();

    // stacktimer::info("wss_internal       %zu", wss_internal);
    // stacktimer::info("wss_persistent     %zu", wss_persistent);
    // stacktimer::info("wss_tmp_preprocess %zu", wss_tmp_preprocess);
    // stacktimer::info("wss_tmp_perform    %zu", wss_tmp_perform);

    stacktimer::pop();

    called_setup = true;
}



template<typename T, typename I>
size_t submatrix_dcsx_ddny<T,I>::get_wss_internal()
{
    if(!called_setup) eslog::error("setup was not called\n");

    return wss_internal;
}



template<typename T, typename I>
size_t submatrix_dcsx_ddny<T,I>::get_wss_persistent()
{
    if(!called_setup) eslog::error("setup was not called\n");

    return wss_persistent;
}



template<typename T, typename I>
size_t submatrix_dcsx_ddny<T,I>::get_wss_tmp_preprocess()
{
    if(!called_setup) eslog::error("setup was not called\n");

    return wss_tmp_preprocess;
}



template<typename T, typename I>
size_t submatrix_dcsx_ddny<T,I>::get_wss_tmp_perform()
{
    if(!called_setup) eslog::error("setup was not called\n");

    return wss_tmp_perform;
}



template<typename T, typename I>
void submatrix_dcsx_ddny<T,I>::set_ws_persistent(void * ws_persistent_)
{
    if(ws_persistent_ == nullptr && wss_persistent > 0) eslog::error("persistent workspace is null\n");
    if(ws_persistent != nullptr) eslog::error("cannot re-set persistent workspace\n");

    ws_persistent = ws_persistent_;
}



template<typename T, typename I>
void submatrix_dcsx_ddny<T,I>::preprocess_submit(void * ws_tmp)
{
    stacktimer::push("submatrix_dcsx_ddny::preprocess_submit");

    if(!called_setup) eslog::error("setup has not been called\n");
    if(called_preprocess) eslog::error("preprocess has already been called\n");
    if(ws_tmp == nullptr && wss_tmp_preprocess > 0) eslog::error("temporary workspace is null\n");

    this->internal_preprocess(ws_tmp);

    stacktimer::pop();

    called_preprocess = true;
}



template<typename T, typename I>
void submatrix_dcsx_ddny<T,I>::perform_submit(void * ws_tmp)
{
    stacktimer::push("submatrix_dcsx_ddny::perform_submit");

    if(!called_preprocess) eslog::error("preprocess has not been called\n");
    if(ws_tmp == nullptr && wss_tmp_perform > 0) eslog::error("temporary workspace is null\n");

    this->internal_perform(ws_tmp);

    stacktimer::pop();
}



#define INSTANTIATE_T_I(T,I) \
template class submatrix_dcsx_ddny<T,I>;

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
