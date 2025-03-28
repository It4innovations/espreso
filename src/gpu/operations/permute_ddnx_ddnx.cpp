
#include "gpu/operations/permute_ddnx_ddnx.h"

#include "basis/utilities/stacktimer.h"
#include "wrappers/cuda/operations/w_cuda_permute_ddnx_ddnx.h"



namespace espreso {
namespace gpu {
namespace operations {



template<typename T, typename I>
std::unique_ptr<permute_ddnx_ddnx<T,I>> permute_ddnx_ddnx<T,I>::make()
{
    #ifdef ESPRESO_USE_WRAPPER_GPU_CUDA
        return std::make_unique<w_cuda_permute_ddnx_ddnx<T,I>>();
    #endif
    eslog::error("wrapper for permute_ddnx_ddnx not available\n");
}



template<typename T, typename I>
void permute_ddnx_ddnx<T,I>::set_handles(gpu::mgm::queue q_)
{
    if(called_set_handles) eslog::error("handles are already set\n");

    q = q_;

    called_set_handles = true;
}



template<typename T, typename I>
void permute_ddnx_ddnx<T,I>::set_matrix_src(MatrixDenseView_new<T> * M_src_)
{
    if(M_src != nullptr) eslog::error("matrix M_src is already set\n");
    if(M_src_ == nullptr) eslog::error("M_src cannot be nullptr\n");

    M_src = M_src_;
}



template<typename T, typename I>
void permute_ddnx_ddnx<T,I>::set_matrix_dst(MatrixDenseView_new<T> * M_dst_)
{
    if(M_dst != nullptr) eslog::error("matrix M_dst is already set\n");
    if(M_dst_ == nullptr) eslog::error("M_dst cannot be nullptr\n");

    M_dst = M_dst_;
}



template<typename T, typename I>
void permute_ddnx_ddnx<T,I>::set_perm_rows(PermutationView_new<I> * perm_rows_)
{
    if(perm_rows != nullptr) eslog::error("perm_rows is already set\n");
    if(perm_rows_ == nullptr) eslog::error("perm_rows cannot be nullptr\n");

    perm_rows = perm_rows_;
}



template<typename T, typename I>
void permute_ddnx_ddnx<T,I>::set_perm_cols(PermutationView_new<I> * perm_cols_)
{
    if(perm_cols != nullptr) eslog::error("perm_cols is already set\n");
    if(perm_cols_ == nullptr) eslog::error("perm_cols cannot be nullptr\n");

    perm_cols = perm_cols_;
}



template<typename T, typename I>
void permute_ddnx_ddnx<T,I>::setup()
{
    stacktimer::push("permute_ddnx_ddnx::setup");

    if(!called_set_handles) eslog::error("handles are not set\n");
    if(M_src == nullptr) eslog::error("source matrix is not sen\n");
    if(M_dst == nullptr) eslog::error("destination matrix is not set\n");
    if(M_src->nrows != M_dst->nrows || M_src->ncols != M_dst->ncols) eslog::error("matrix sizes dont match\n");
    if(M_src->order != M_dst->order) eslog::error("matrix orders must match\n");
    if(perm_rows != nullptr && perm_rows->size != M_src->nrows) eslog::error("wrong row perm size\n");
    if(perm_cols != nullptr && perm_cols->size != M_src->ncols) eslog::error("wrong col perm size\n");

    perm_primary = perm_rows;
    perm_secdary = perm_cols;
    if(M_src->order == 'C') {
        std::swap(perm_primary, perm_secdary);
    }

    this->internal_setup();

    // stacktimer::info("wss_tmp_perform    %zu", wss_tmp_perform);

    stacktimer::pop();

    called_setup = true;
}



template<typename T, typename I>
size_t permute_ddnx_ddnx<T,I>::get_wss_tmp_perform()
{
    if(!called_setup) eslog::error("setup was not called\n");

    return wss_tmp_perform;
}



template<typename T, typename I>
void permute_ddnx_ddnx<T,I>::perform_submit(void * ws_tmp)
{
    stacktimer::push("permute_ddnx_ddnx::perform_submit");

    if(!called_setup) eslog::error("setup was not called\n");
    if(ws_tmp == nullptr && wss_tmp_perform > 0) eslog::error("temporary workspace is null\n");

    this->internal_perform(ws_tmp);

    stacktimer::pop();
}



#define INSTANTIATE_T_I(T,I) \
template class permute_ddnx_ddnx<T,I>;

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

