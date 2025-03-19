
#include "gpu/operations/auxliary/gpu_trsm_trirhs_chunk_splitfactor.h"



namespace espreso {
namespace math {
namespace operations {



template<typename T, typename I>
void gpu_trsm_trirhs_chunk_splitfactor<T,I>::set_config(config cfg_)
{
    if(called_set_config) eslog::error("config has already been set\n");

    cfg = cfg_;

    called_set_config = true;
}



template<typename T, typename I>
void gpu_trsm_trirhs_chunk_splitfactor<T,I>::set_range(size_t k_start_, size_t k_end_)
{
    if(called_set_range) eslog::error("range has already been set\n");

    k_start = k_start_;
    k_end = k_end_;
    k_size = k_end - k_start;

    called_set_range = true;
}



template<typename T, typename I>
void gpu_trsm_trirhs_chunk_splitfactor<T,I>::set_handles(gpu::mgm::queue q_, gpu::spblas::handle handle_spblas_, gpu::dnblas::handle handle_dnblas_)
{
    if(called_set_handles) eslog::error("handles have already been set\n");

    q = q_;
    handle_spblas = handle_spblas_;
    handle_dnblas = handle_dnblas_;

    called_set_handles = true;
}



template<typename T, typename I>
void gpu_trsm_trirhs_chunk_splitfactor<T,I>::set_matrix_h_L(MatrixCsxView_new<T,I> * h_L_)
{
    if(h_L != nullptr) eslog::error("matrix h_L is already set\n");
    if(h_L_ == nullptr) eslog::error("h_L cannot be nullptr\n");

    h_L = h_L_;
}



template<typename T, typename I>
void gpu_trsm_trirhs_chunk_splitfactor<T,I>::set_matrix_d_X(MatrixDenseView_new<T> d_X_)
{
    if(d_X != nullptr) eslog::error("matrix d_X is already set\n");
    if(d_X_ == nullptr) eslog::error("d_X cannot be nullptr\n");

    d_X = d_X_;
}



template<typename T, typename I>
void gpu_trsm_trirhs_chunk_splitfactor<T,I>::set_h_X_rowtrails(VectorDenseView_new<I> h_X_rowtrails_)
{
    if(h_X_rowtrails != nullptr) eslog::error("X rowtrails are already set\n");
    if(h_X_rowtrails_ == nullptr) eslog::error("X rowtrails cannot be nullptr\n");

    h_X_rowtrails = h_X_rowtrails_;
}



template<typename T, typename I>
void gpu_trsm_trirhs_chunk_splitfactor<T,I>::setup()
{
    if(!called_set_config) eslog::error("config is not set\n");
    if(!called_set_range) eslog::error("range is not set\n");
    if(!called_set_handles) eslog::error("handles are not set\n");
    if(called_setup) eslog::error("setup has already been called\n");
    if(h_L == nullptr) eslog::error("matrix L is not set\n");
    if(d_X == nullptr) eslog::error("matrix X is not set\n");
    if(h_X_rowtrails == nullptr) eslog::error("X rowtrails are not set\n");
    if(L->nrows != L->ncols) eslog::error("L is not square\n");
    if(L->nrows != X->nrows) eslog::error("incompatible matrices\n");
    if(h_X_rowtrails->size != X->nrows) eslog::error("wrong X rowtrails size\n");
    if(h_L->prop.uplo != 'L') eslog::error("matrix L must have uplo=L\n");

    ator_ws_persistent = std::make_unique<AllocatorArena_new>(false, true, gpu::mgm::get_natural_pitch_align());
    ator_ws_tmp_linear = std::make_unique<AllocatorArena_new>(false, true, gpu::mgm::get_natural_pitch_align());
    ator_ws_tmp_overlap = std::make_unique<AllocatorSinglePointer_new>(false, true, gpu::mgm::get_natural_pitch_align());

    op_h_submatrix_L_top.set_matrix_src(h_L);
    op_h_submatrix_L_top.set_bounds(k_start, k_end, k_start, k_end);
    op_h_submatrix_L_top.setup();
    size_t nnz_L_top = op_h_submatrix_L_top.get_output_matrix_nnz();
    h_sub_L_top_sp.set(k_size, k_size, nnz_L_top, cfg.trsm_factor_order, AllocatorHostPinend_new::get_singleton());
    h_sub_L_top_sp.prop.uplo = h_L->prop.uplo;
    h_sub_L_top_sp.prop.diag = h_L->prop.diag;
    op_h_submatrix_L_top.set_matrix_dst(h_sub_L_top_sp);

    op_h_submatrix_L_bot.set_matrix_src(h_L);
    op_h_submatrix_L_bot.set_bounds(k_end, h_L->nrows, k_start, k_end);
    op_h_submatrix_L_bot.setup();
    size_t nnz_L_bot = op_h_submatrix_L_bot.get_output_matrix_nnz();
    h_sub_L_bot_sp.set(h_L.nrows - k_end, k_size, nnz_L_bot, cfg.gemm_factor_order, AllocatorHostPinned_new::get_singleton());
    op_h_submatrix_L_bot.set_matrix_dst(h_sub_L_bot_sp);

    rhs_end = h_X_rowtrails.vals[k_end - 1] + 1;

    d_sub_X_top.set_view(k_size, rhs_end, d_X->order, nullptr);
    op_submatrix_d_X_top.set_matrix_src(d_X);
    op_submatrix_d_X_top.set_matrix_dst(&d_sub_X_top);
    op_submatrix_d_X_top.set_bounds(k_start, k_end, 0, rhs_end);

    d_sub_X_bot.set_view(d_X->nrows - k_end, rhs_end, d_X->order, nullptr);
    op_submatrix_d_X_bot.set_matrix_src(d_X);
    op_submatrix_d_X_bot.set_matrix_dst(&d_sub_X_bot);
    op_submatrix_d_X_bot.set_bounds(k_end, X->nrows, 0, rhs_end);

    op_trsm.set_config(cfg.trsm_factor_spdn);
    op_trsm.set_handles(q, handle_spblas, handle_dnblas);
    op_trsm.set_matrix_h_A(&h_sub_L_top_sp);
    op_trsm.set_matrix_d_X(&d_sub_X_top);
    op_trsm.set_matrix_d_B(&d_sub_X_top);
    op_trsm.setup();
    wss_internal += op_trsm.get_wss_internal();
    wss_persistent += op_trsm.get_wss_persistent();
    wss_tmp_preprocess_overlap = std::max(wss_tmp_preprocess_overlap, op_trsm.get_wss_tmp_preprocess());
    wss_tmp_perform_overlap = std::max(wss_tmp_perform_overlap, op_trsm.get_wss_tmp_perform());

    if(d_sub_X_bot.nrows > 0) {
        bool prune_rows = ((cfg.gemm_factor_prune == 'R') || (cfg.gemm_factor_prune == 'A'));
        bool prune_cols = ((cfg.gemm_factor_prune == 'C') || (cfg.gemm_factor_prune == 'A'));
        op_gemm.set_config(cfg.gemm_factor_spdn, prune_rows, prune_cols);
        op_gemm.set_handles(q, handle_spblas, handle_dnblas);
        op_gemm.set_matrix_h_A(&h_sub_L_bot_sp);
        op_gemm.set_matrix_d_B(d_sub_X_top);
        op_gemm.set_matrix_d_C(d_sub_X_bot);
        op_gemm.set_coefficients(T{-1}, T{1});
        op_gemm.setup();
        wss_internal += op_gemm.get_wss_internal();
        wss_persistent += op_gemm.get_wss_persistent();
        wss_tmp_preprocess_overlap = std::max(wss_tmp_preprocess_overlap, op_gemm.get_wss_tmp_preprocess());
        wss_tmp_perform_overlap = std::max(wss_tmp_perform_overlap, op_gemm.get_wss_tmp_perform());
    }

    wss_tmp_preprocess_linear = ((wss_tmp_preprocess_linear - 1) / ator_ws_tmp_linear->get_align() + 1) * ator_ws_tmp_linear->get_align();
    wss_tmp_perform_linear = ((wss_tmp_perform_linear - 1) / ator_ws_tmp_linear->get_align() + 1) * ator_ws_tmp_linear->get_align();

    wss_tmp_preprocess = wss_tmp_preprocess_linear + wss_tmp_preprocess_overlap;
    wss_tmp_perform = wss_tmp_perform_linear + wss_tmp_perform_overlap;

    called_setup = true;
}



template<typename T, typename I>
size_t gpu_trsm_trirhs_chunk_splitfactor<T,I>::get_wss_internal()
{
    if(!called_setup) eslog::error("setup was not called\n");

    return wss_internal;
}



template<typename T, typename I>
size_t gpu_trsm_trirhs_chunk_splitfactor<T,I>::get_wss_persistent()
{
    if(!called_setup) eslog::error("setup was not called\n");

    return wss_persistent;
}



template<typename T, typename I>
size_t gpu_trsm_trirhs_chunk_splitfactor<T,I>::get_wss_tmp_preprocess()
{
    if(!called_setup) eslog::error("setup was not called\n");

    return wss_tmp_preprocess;
}



template<typename T, typename I>
size_t gpu_trsm_trirhs_chunk_splitfactor<T,I>::get_wss_tmp_perform()
{
    if(!called_setup) eslog::error("setup was not called\n");

    return wss_tmp_perform;
}



template<typename T, typename I>
void gpu_trsm_trirhs_chunk_splitfactor<T,I>::set_ws_persistent(void * ws_persistent_)
{
    if(ws_persistent_ == nullptr && wss_persistent > 0) eslog::error("persistent workspace is null\n");
    if(ws_persistent != nullptr) eslog::error("cannot re-set persistent workspace\n");

    ws_persistent = ws_persistent_;
}



template<typename T, typename I>
void gpu_trsm_trirhs_chunk_splitfactor<T,I>::preprocess_submit(void * ws_tmp)
{
    if(!called_setup) eslog::error("setup has not been called\n");
    if(called_preprocess) eslog::error("preprocess has already been called\n");
    if(ws_tmp == nullptr && wss_tmp_preprocess > 0) eslog::error("temporary workspace is null\n");

    ator_ws_persistent.set(ws_persistent, wss_persistent);

    ator_ws_tmp_linear.set(ws_tmp, wss_tmp_preprocess_linear);
    ator_ws_tmp_overlap.set(ws_tmp + wss_tmp_preprocess_linear, wss_tmp_preprocess_overlap);

    h_sub_L_top_sp.alloc();
    op_h_submatrix_L_top.set_matrix_dst(&h_sub_L_top_sp);
    op_h_submatrix_L_top.perform_pattern();

    h_sub_L_bot_sp.alloc();
    op_h_submatrix_L_bot.set_matrix_dst(&h_sub_L_bot_sp);
    op_h_submatrix_L_bot.perform_pattern();

    op_trsm.set_ws_persistent(ator_ws_persistent->alloc(op_trsm.get_wss_persistent()));
    op_trsm.preprocess_submit(ator_ws_tmp_overlap->alloc(op_trsm.get_wss_tmp_preprocess()));

    if(d_sub_X_bot.nrows > 0) {
        op_gemm.set_ws_persistent(ator_ws_persistent->alloc(op_gemm.get_wss_persistent()));
        op_gemm.preprocess_submit(ator_ws_tmp_overlap->alloc(op_trsm.get_wss_tmp_preprocess()));
    }

    ator_ws_tmp_linear.unset();
    ator_ws_tmp_overlap.unset();

    called_preprocess = true;
}



template<typename T, typename I>
void gpu_trsm_trirhs_chunk_splitfactor<T,I>::perform_submit(void * ws_tmp)
{
    if(!called_preprocess) eslog::error("preprocess has not been called\n");
    if(ws_tmp == nullptr && wss_tmp_perform > 0) eslog::error("temporary workspace is null\n");

    ator_ws_tmp_linear.set(ws_tmp, wss_tmp_perform_linear);
    ator_ws_tmp_overlap.set(ws_tmp + wss_tmp_perform_linear, wss_tmp_perform_overlap);

    op_h_submatrix_L_top.perform_values();

    op_h_submatrix_L_bot.perform_values();

    op_submatrix_d_X_top.perform();

    op_submatrix_d_X_bot.perform();

    op_trsm.perform_submit(ator_ws_tmp_overlap->alloc(op_trsm.get_wss_tmp_perform()));

    if(d_sub_X_bot.nrows > 0) {
        op_gemm.perform_submit(ator_ws_tmp_overlap->alloc(op_trsm.get_wss_tmp_perform()));
    }

    ator_ws_tmp_linear.unset();
    ator_ws_tmp_overlap.unset();
}



#define INSTANTIATE_T_I(T,I) \
template class gpu_trsm_trirhs_chunk_splitfactor<T,I>;

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
