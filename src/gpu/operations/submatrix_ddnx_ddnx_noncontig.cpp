
#include "gpu/operations/submatrix_ddnx_ddnx_noncontig.h"

#include "basis/utilities/stacktimer.h"
#include "wrappers/cuda/operations/w_cuda_submatrix_ddnx_ddnx_noncontig.h"
#include "wrappers/rocm/operations/w_rocm_submatrix_ddnx_ddnx_noncontig.h"
#include "wrappers/oneapi/operations/w_oneapi_submatrix_ddnx_ddnx_noncontig.h"



namespace espreso {
namespace gpu {
namespace operations {



template<typename T, typename I>
std::unique_ptr<submatrix_ddnx_ddnx_noncontig<T,I>> submatrix_ddnx_ddnx_noncontig<T,I>::make()
{
    #ifdef ESPRESO_USE_WRAPPER_GPU_CUDA
        return std::make_unique<w_cuda_submatrix_ddnx_ddnx_noncontig<T,I>>();
    #endif
    #ifdef ESPRESO_USE_WRAPPER_GPU_ROCM
        return std::make_unique<w_rocm_submatrix_ddnx_ddnx_noncontig<T,I>>();
    #endif
    #ifdef ESPRESO_USE_WRAPPER_GPU_ONEAPI
        return std::make_unique<w_oneapi_submatrix_ddnx_ddnx_noncontig<T,I>>();
    #endif
    eslog::error("wrapper for submatrix_ddnx_ddnx_noncontig not available\n");
}



template<typename T, typename I>
void submatrix_ddnx_ddnx_noncontig<T,I>::set_handles(gpu::mgm::queue q_)
{
    if(called_set_handles) eslog::error("handles are already set\n");

    q = q_;

    called_set_handles = true;
}



template<typename T, typename I>
void submatrix_ddnx_ddnx_noncontig<T,I>::set_matrix_src(MatrixDenseView_new<T> * d_M_src_)
{
    if(d_M_src != nullptr) eslog::error("source matrix is already set\n");
    if(d_M_src_ == nullptr) eslog::error("source matrix cannot be nullptr\n");

    d_M_src = d_M_src_;
}



template<typename T, typename I>
void submatrix_ddnx_ddnx_noncontig<T,I>::set_matrix_dst(MatrixDenseView_new<T> * d_M_dst_)
{
    if(d_M_dst != nullptr) eslog::error("destination matrix is already set\n");
    if(d_M_dst_ == nullptr) eslog::error("destination matrix cannot be nullptr\n");

    d_M_dst = d_M_dst_;
}



template<typename T, typename I>
void submatrix_ddnx_ddnx_noncontig<T,I>::set_row_map(VectorDenseView_new<I> * d_row_map_)
{
    if(called_set_row_map) eslog::error("row map is already set\n");

    d_row_map = d_row_map_;

    called_set_row_map = true;
}



template<typename T, typename I>
void submatrix_ddnx_ddnx_noncontig<T,I>::set_col_map(VectorDenseView_new<I> * d_col_map_)
{
    if(called_set_col_map) eslog::error("col map is already set\n");

    d_col_map = d_col_map_;

    called_set_col_map = true;
}



template<typename T, typename I>
void submatrix_ddnx_ddnx_noncontig<T,I>::perform_submit()
{
    stacktimer::push("submatrix_ddnx_ddnx_noncontig::perform_submit");

    if(!called_set_handles) eslog::error("handles are not set\n");
    if(d_M_src == nullptr) eslog::error("source matrix is not set\n");
    if(d_M_dst == nullptr) eslog::error("destination matrix is not set\n");
    if(!d_M_src->ator->is_data_accessible_gpu()) eslog::error("source matrix must be gpu-accessible\n");
    if(!d_M_dst->ator->is_data_accessible_gpu()) eslog::error("destination matrix must be gpu-accessible\n");
    if(d_row_map != nullptr && !d_row_map->ator->is_data_accessible_gpu()) eslog::error("d_row_map must be gpu-accessible\n");
    if(d_col_map != nullptr && !d_col_map->ator->is_data_accessible_gpu()) eslog::error("d_col_map must be gpu-accessible\n");
    if(d_M_src->order != d_M_dst->order) eslog::error("matrix orders do no tmatch");
    if(d_row_map != nullptr && d_row_map->size != d_M_dst->nrows) eslog::error("incompatible row sizes\n");
    if(d_col_map != nullptr && d_col_map->size != d_M_dst->ncols) eslog::error("incompatible col sizes\n");

    this->internal_perform();

    stacktimer::pop();
}



template<typename T, typename I>
void submatrix_ddnx_ddnx_noncontig<T,I>::submit_all(gpu::mgm::queue q, MatrixDenseView_new<T> * d_M_src, MatrixDenseView_new<T> * d_M_dst, VectorDenseView_new<I> * d_row_map, VectorDenseView_new<I> * d_col_map)
{
    auto instance = submatrix_ddnx_ddnx_noncontig<T,I>::make();
    instance->set_handles(q);
    instance->set_matrix_src(d_M_src);
    instance->set_matrix_dst(d_M_dst);
    instance->set_row_map(d_row_map);
    instance->set_col_map(d_col_map);
    instance->perform_submit();
}



#define INSTANTIATE_T_I(T,I) \
template class submatrix_ddnx_ddnx_noncontig<T,I>;

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
