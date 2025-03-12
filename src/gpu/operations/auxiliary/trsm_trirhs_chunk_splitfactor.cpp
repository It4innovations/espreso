
#include "gpu/operations/auxliary/trsm_trirhs_chunk_splitfactor.h"



namespace espreso {
namespace math {
namespace operations {



template<typename T, typename I>
void trsm_trirhs_chunk_splitfactor<T,I>::set_config(config cfg_)
{
    if(called_set_config) eslog::error("config has already been set\n");

    cfg = cfg_;

    called_set_config = true;
}



template<typename T, typename I>
void trsm_trirhs_chunk_splitfactor<T,I>::set_range(size_t k_start_, size_t k_end_)
{
    if(called_set_range) eslog::error("range has already been set\n");

    k_start = k_start_;
    k_end = k_end_;
    k_size = k_end - k_start;

    called_set_range = true;
}



template<typename T, typename I>
void trsm_trirhs_chunk_splitfactor<T,I>::set_handles(gpu::mgm::queue q_, gpu::spblas::handle spblas_handle_, gpu::dnblas::handle dnblas_handle_)
{
    if(called_set_handles) eslog::error("handles have already been set\n");

    q = q_;
    handle_spblas = spblas_handle_;
    handle_dnblas = dnblas_handle_;

    called_set_handles = true;
}



template<typename T, typename I>
void trsm_trirhs_chunk_splitfactor<T,I>::set_matrix_h_L(MatrixCsxView_new<T,I> h_L_)
{
    if(called_set_L && !MatrixCsxView_new<T,I>::are_interchangable(h_L, h_L_)) eslog::error("invalid replacement for matrix L\n");

    h_L = h_L_;

    op_h_submatrix_L_top.set_matrix_src(&h_L);
    op_h_submatrix_L_bot.set_matrix_src(&h_L);

    called_set_L = true;
}



template<typename T, typename I>
void trsm_trirhs_chunk_splitfactor<T,I>::set_matrix_d_X(MatrixDenseView_new<T> d_X_)
{
    if(called_set_X && !MatrixCsxView_new<T,I>::are_interchangable(d_X, d_X_)) eslog::error("invalid replacement for matrix X\n");

    d_X = d_X_;

    op_submatrix_d_X_top.set_matrix_src(&d_X);
    op_submatrix_d_X_bot.set_matrix_src(&d_X);

    called_set_X = true;
}



template<typename T, typename I>
void trsm_trirhs_chunk_splitfactor<T,I>::set_h_X_rowtrails(VectorDenseView_new<I> h_X_rowtrails_)
{
    if(called_set_X_rowtrails) eslog::error("X rowtrails have already been set\n");

    h_X_rowtrails = h_X_rowtrails_;

    called_set_X_rowtrails = true;
}



template<typename T, typename I>
void trsm_trirhs_chunk_splitfactor<T,I>::setup()
{
    if(!called_set_config) eslog::error("config is not set\n");
    if(!called_set_range) eslog::error("range is not set\n");
    if(!called_set_handles) eslog::error("handles are not set\n");
    if(!called_set_L) eslog::error("matrix L is not set\n");
    if(!called_set_X) eslog::error("matrix X is not set\n");
    if(!called_set_X_rowtrails) eslog::error("X rowtrails are not set\n");
    if(L.nrows != L.ncols) eslog::error("L is not square\n");
    if(L.nrows != X.nrows) eslog::error("incompatible matrices\n");
    if(h_X_rowtrails.size != X.nrows) eslog::error("wrong X rowtrails size\n");

    op_h_submatrix_L_top.set_bounds(k_start, k_end, k_start, k_end);
    op_h_submatrix_L_top.setup();
    size_t nnz_L_top = op_h_submatrix_L_top.get_output_matrix_nnz();

    op_h_submatrix_L_bot.set_bounds(k_end, h_L.nrows, k_start, k_end);
    op_h_submatrix_L_bot.setup();
    size_t nnz_L_bot = op_h_submatrix_L_bot.get_output_matrix_nnz();

    h_sub_L_top_sp.set(k_size, k_size, nnz_L_top, cfg.trsm_factor_order, AllocatorCPU_new::get_singleton());
    h_sub_L_bot_sp.set(h_L.nrows - k_end, k_size, nnz_L_bot, cfg.gemm_factor_order, AllocatorCPU_new::get_singleton());

    d_sub_L_top_sp.set(h_sub_L_top_sp.nrows, h_sub_L_top_sp.ncols, h_sub_L_top_sp.nnz, h_sub_L_top_sp.order, ator_ws_persistent.get());
    d_sub_L_bot_sp.set(h_sub_L_bot_sp.nrows, h_sub_L_bot_sp.ncols, h_sub_L_bot_sp.nnz, h_sub_L_bot_sp.order, ator_ws_persistent.get());



    // todo

    called_setup = true;
}



template<typename T, typename I>
size_t trsm_trirhs_chunk_splitfactor<T,I>::get_wss_internal()
{
    if(!called_setup) eslog::error("setup was not called\n");

    return wss_internal;
}



template<typename T, typename I>
size_t trsm_trirhs_chunk_splitfactor<T,I>::get_wss_persistent()
{
    if(!called_setup) eslog::error("setup was not called\n");

    return wss_persistent;
}



template<typename T, typename I>
size_t trsm_trirhs_chunk_splitfactor<T,I>::get_wss_tmp_preprocess()
{
    if(!called_setup) eslog::error("setup was not called\n");

    return wss_tmp_preprocess;
}



template<typename T, typename I>
size_t trsm_trirhs_chunk_splitfactor<T,I>::get_wss_tmp_perform()
{
    if(!called_setup) eslog::error("setup was not called\n");

    return wss_tmp_perform;
}



template<typename T, typename I>
void trsm_trirhs_chunk_splitfactor<T,I>::set_ws_persistent(void * ws_persistent_)
{
    if(ws_persistent_ == nullptr && wss_persistent > 0) eslog::error("persistent workspace is null\n");
    if(ws_persistent != nullptr) eslog::error("cannot re-set persistent workspace\n");

    ws_persistent = ws_persistent_;
}



template<typename T, typename I>
void trsm_trirhs_chunk_splitfactor<T,I>::preprocess_submit(void * ws_tmp)
{
    if(!called_setup) eslog::error("setup has not been called\n");
    if(called_preprocess) eslog::error("preprocess has already been called\n");
    if(ws_tmp == nullptr && wss_tmp_preprocess > 0) eslog::error("temporary workspace is null\n");

    ator_ws_persistent = std::make_unique<AllocatorArena_new>(ws_persistent, wss_persistent, gpu::mgm::get_natural_pitch_align(), false, true);

    h_sub_L_top_sp.alloc();
    op_h_submatrix_L_top.set_matrix_dst(&h_sub_L_top_sp);
    op_h_submatrix_L_top.perform_pattern();
    d_sub_L_top_sp.alloc();
    gpu::mgm::copy_submit(q, h_sub_L_top_sp, d_sub_L_top_sp, true, false);

    h_sub_L_bot_sp.alloc();
    op_h_submatrix_L_bot.set_matrix_dst(&h_sub_L_bot_sp);
    op_h_submatrix_L_bot.perform_pattern();
    d_sub_L_bot_sp.alloc();
    gpu::mgm::copy_submit(q, h_sub_L_bot_sp, d_sub_L_bot_sp, true, false);


    



    size_t rhs_end = h_X_rowtrails.vals[k_end - 1] + 1;

    d_sub_X_top.set_view(k_size, rhs_end, d_X.order, nullptr);

    op_submatrix_d_X_top.set_matrix_src(&d_X);
    op_submatrix_d_X_top.set_matrix_dst(&d_sub_X_top);
    op_submatrix_d_X_top.set_bounds(k_start, k_end, 0, rhs_end);

    op_d_trsm_sp.set_handles(q, handle_spblas);
    op_d_trsm_sp.set_matrix_A(d_sub_L_top_sp);
    op_d_trsm_sp.set_matrix_X(d_sub_X_top);
    op_d_trsm_sp.set_matrix_B(d_sub_X_top);
    op_d_trsm_sp.setup();

    // todo
}



template<typename T, typename I>
void trsm_trirhs_chunk_splitfactor<T,I>::update_submit()
{
    if(!called_preprocess) eslog::error("preprocess has not been called\n");

    // todo
}



template<typename T, typename I>
void trsm_trirhs_chunk_splitfactor<T,I>::perform_submit(void * ws_tmp)
{
    if(!called_preprocess) eslog::error("preprocess has not been called\n");
    if(ws_tmp == nullptr && wss_tmp_perform > 0) eslog::error("temporary workspace is null\n");

    op_submatrix_d_X_top.perform();

    // trsm top (sparse/dense)
    // optionally convert to dense
    // trsm Ltop Xtop

    op_h_submatrix_L_top.perform_values();
    gpu::mgm::copy_submit(q, h_sub_L_top_sp, d_sub_L_top_sp, false, true);

    op_h_submatrix_L_bot.perform_values();
    gpu::mgm::copy_submit(q, h_sub_L_bot_sp, d_sub_L_bot_sp, false, true);

    if(cfg.trsm_factor_spdn == 'S') {

    }
    if(cfg.trsm_factor_spdn == 'D') {
        
    }


    // gemm bot (sparse/dense, prune/noprune)
}



#define INSTANTIATE_T_I(T,I) \
template class trsm_trirhs_chunk_splitfactor<T,I>;

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
