
#include "gpu/operations/gemm_hcsx_ddny_ddnz_prune.h"



namespace espreso {
namespace gpu {
namespace operations {



template<typename T, typename I>
void gemm_hcsx_ddny_ddnz_prune<T,I>::set_config(char spdn_A_, bool prune_rows_, bool prune_cols_)
{
    if(called_set_config) eslog::error("config are already set\n");

    spdn_A = spdn_A_;
    prune_rows = prune_rows_;
    prune_cols = prune_cols_;

    called_set_config = true;
}



template<typename T, typename I>
void gemm_hcsx_ddny_ddnz_prune<T,I>::set_handles(gpu::mgm::queue q_, gpu::spblas::handle handle_spblas, gpu::dnblas::handle dnblas_handle_)
{
    if(called_set_handles) eslog::error("handles are already set\n");

    q = q_;
    handle_spblas = handle_spblas_;
    handle_dnblas = handle_dnblas_;

    called_set_handles = true;
}



template<typename T, typename I>
void gemm_hcsx_ddny_ddnz_prune<T,I>::set_matrix_h_A(MatrixCsxView_new<T,I> h_A_)
{
    if(called_set_A && !MatrixCsxView_new<T,I>::are_interchangable(h_A, h_A_)) eslog::error("invalid replacement for matrix A\n");

    h_A = h_A_;

    op_h_prune_A.set_matrix_src(&h_A);

    called_set_A = true;
}



template<typename T, typename I>
void gemm_hcsx_ddny_ddnz_prune<T,I>::set_matrix_d_B(MatrixDenseView_new<T> d_B_)
{
    if(called_set_B && !MatrixDenseView_new<T>::are_interchangable(d_B, d_B_)) eslog::error("invalid replacement for matrix B\n");

    d_B = d_B_;

    called_set_B = true;
}



template<typename T, typename I>
void gemm_hcsx_ddny_ddnz_prune<T,I>::set_matrix_d_C(MatrixDenseView_new<T> d_C_)
{
    if(called_set_C && !MatrixDenseView_new<T>::are_interchangable(d_C, d_C_)) eslog::error("invalid replacement for matrix C\n");

    d_C = d_C_;

    called_set_C = true;
}



template<typename T, typename I>
void gemm_hcsx_ddny_ddnz_prune<T,I>::set_coefficients(T alpha_, T beta_)
{
    alpha = alpha_;
    beta = beta_;
}



template<typename T, typename I>
void gemm_hcsx_ddny_ddnz_prune<T,I>::setup()
{
    if(!called_set_config) eslog::error("config is not set\n");
    if(!called_set_handles) eslog::error("handles are not set\n");
    if(!called_set_A) eslog::error("matrix A is not set\n");
    if(!called_set_B) eslog::error("matrix B is not set\n");
    if(!called_set_C) eslog::error("matrix C is not set\n");
    if(called_setup) eslog::error("setup has already been called\n");
    if(A.nrows != C.nrows || B.ncols != C.ncols || A.ncols != B.nrows) eslog::error("incompatible matrices\n");

    if(prune_rows || prune_cols) {
        op_h_prune_A.set_pruning_mode(prune_rows, prune_cols);
        op_h_prune_A.setup();

        m = op_h_prune_A.get_dst_matrix_nrows();
        n = B->ncols;
        k = op_h_prune_A.get_dst_matrix_ncols();

        if(prune_rows) {
            pruned_rows.set(m, AllocatorCPU_new::get_singleton());
            pruned_rows.alloc();
            op_h_prune_A.set_vector_pruned_rows(&pruned_rows);
        }
        if(prune_cols) {
            pruned_cols.set(k, AllocatorCPU_new::get_singleton());
            pruned_cols.alloc();
            op_h_prune_A.set_vector_pruned_cols(&pruned_cols);
        }

        op_prune_A.preprocess2();
        
        op_prune_A.set_matrix_dst_sp(&h_A_pruned);

        h_A_to_use = &h_A_pruned;
    }
    else {
        h_A_to_use = &h_A;
    }
    
    d_A_pruned_sp.set(h_A_to_use->nrows, h_A_to_use->ncols, h_A_to_use->nnz, h_A_to_use->order, );
    
    
    
    // todo

    called_setup = true;
}



template<typename T, typename I>
size_t gemm_hcsx_ddny_ddnz_prune<T,I>::get_wss_internal()
{
    return wss_internal;
}



template<typename T, typename I>
size_t gemm_hcsx_ddny_ddnz_prune<T,I>::get_wss_persistent()
{
    return wss_persistent;
}



template<typename T, typename I>
size_t gemm_hcsx_ddny_ddnz_prune<T,I>::get_wss_tmp_preprocess()
{
    return wss_tmp_preprocess;
}



template<typename T, typename I>
size_t gemm_hcsx_ddny_ddnz_prune<T,I>::get_wss_tmp_perform()
{
    return wss_tmp_perform;
}



template<typename T, typename I>
void gemm_hcsx_ddny_ddnz_prune<T,I>::set_ws_persistent(void * ws_persistent_)
{
    if(ws_persistent_ == nullptr && wss_persistent > 0) eslog::error("persistent workspace is null\n");
    if(ws_persistent != nullptr) eslog::error("cannot re-set persistent workspace\n");

    ws_persistent = ws_persistent_;
}



template<typename T, typename I>
void gemm_hcsx_ddny_ddnz_prune<T,I>::preprocess_submit(void * ws_tmp)
{
    if(!called_setup) eslog::error("setup has not been called\n");
    if(called_preprocess) eslog::error("preprocess has already been called\n");
    if(ws_tmp == nullptr && wss_tmp_preprocess > 0) eslog::error("temporary workspace is null\n");

    // todo

    called_preprocess = true;
}



template<typename T, typename I>
void gemm_hcsx_ddny_ddnz_prune<T,I>::perform_submit(void * ws_tmp)
{
    if(!called_preprocess) eslog::error("preprocess has not been called\n");
    if(ws_tmp == nullptr && wss_tmp_perform > 0) eslog::error("temporary workspace is null\n");

    // todo
}



#define INSTANTIATE_T_I(T,I) \
template class gemm_hcsx_ddny_ddnz_prune<T,I>;

    #define INSTANTIATE_T(T) \
    INSTANTIATE_T_I(T,int32_t) \
    /* INSTANTIATE_T_I(T,int64_t) */

        #define INSTANTIATE \
        /* INSTANTIATE_T(float) */ \
        INSTANTIATE_T(double) \
        /* INSTANTIATE_T(std::complex<float>) */ \
        /* INSTANTIATE_T(std::complex<double>) */

            INSTANTIATE

        #undef INSTANTIATE
    #undef INSTANTIATE_T
#undef INSTANTIATE_T_I



}
}
}
