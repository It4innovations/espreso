
#include "gpu/operations/auxiliary/gpu_herk_tri_chunk_stairs.h"

#include "basis/utilities/stacktimer.h"



namespace espreso {
namespace gpu {
namespace operations {



template<typename T, typename I>
void gpu_herk_tri_chunk_stairs<T,I>::set_range(size_t n_start_, size_t n_end_)
{
    if(called_set_range) eslog::error("range is already set\n");

    n_start = n_start_;
    n_end = n_end_;

    called_set_range = true;
}



template<typename T, typename I>
void gpu_herk_tri_chunk_stairs<T,I>::set_handles(gpu::mgm::queue q_, gpu::dnblas::handle handle_dnblas_)
{
    if(called_set_handles) eslog::error("handles are already set\n");

    q = q_;
    handle_dnblas = handle_dnblas_;

    called_set_handles = true;
}



template<typename T, typename I>
void gpu_herk_tri_chunk_stairs<T,I>::set_matrix_d_A_left(MatrixDenseView_new<T> * d_A_left_)
{
    if(d_A_left != nullptr) eslog::error("d_A_left is already set\n");
    if(d_A_left_ == nullptr) eslog::error("d_A_left cannot be nullptr\n");

    d_A_left = d_A_left_;
}



template<typename T, typename I>
void gpu_herk_tri_chunk_stairs<T,I>::set_matrix_d_A_top(MatrixDenseView_new<T> * d_A_top_)
{
    if(d_A_top != nullptr) eslog::error("d_A_top is already set\n");
    if(d_A_top_ == nullptr) eslog::error("d_A_top cannot be nullptr\n");

    d_A_top = d_A_top_;
}



template<typename T, typename I>
void gpu_herk_tri_chunk_stairs<T,I>::set_matrix_d_C(MatrixDenseView_new<T> * d_C_)
{
    if(d_C != nullptr) eslog::error("d_C is already set\n");
    if(d_C_ == nullptr) eslog::error("d_C cannot be nullptr\n");

    d_C = d_C_;
}



template<typename T, typename I>
void gpu_herk_tri_chunk_stairs<T,I>::set_h_A_pivots(VectorDenseView_new<I> * h_A_pivots_)
{
    if(h_A_pivots != nullptr) eslog::error("h_A_pivots is already set\n");
    if(h_A_pivots_ == nullptr) eslog::error("h_A_pivots cannot be nullptr\n");

    h_A_pivots = h_A_pivots_;
}



template<typename T, typename I>
void gpu_herk_tri_chunk_stairs<T,I>::set_coefficients(Treal alpha_, Treal beta_)
{
    alpha = alpha_;
    beta = beta_;
}



template<typename T, typename I>
void gpu_herk_tri_chunk_stairs<T,I>::setup()
{
    stacktimer::push("gpu_herk_tri_chunk_stairs::setup");

    if(utils::is_complex<T>()) eslog::error("complex numbers not supported\n");
    if(!called_set_range) eslog::error("range is not set\n");
    if(!called_set_handles) eslog::error("handles are not set\n");
    if(called_setup) eslog::error("setup was already called\n");
    if(d_A_left == nullptr) eslog::error("d_A_left is not set\n");
    if(d_A_top == nullptr) eslog::error("d_A_top is not set\n");
    if(d_C == nullptr) eslog::error("d_C is not set\n");
    if(h_A_pivots == nullptr) eslog::error("h_A_pivots is not set\n");
    if(d_A_left->nrows != d_A_top->ncols || d_A_left->ncols != d_A_top->nrows) eslog::error("wrong matrices A\n");
    if(d_A_left->order == d_A_top->order) eslog::error("wrong matrices A\n");
    if(d_A_left->vals != d_A_top->vals) eslog::error("wrong matrices A\n");
    if(d_C->nrows != d_C->ncols) eslog::error("d_C must be square\n");
    if(d_C->prop.uplo != 'U' && d_C->prop.uplo != 'L') eslog::error("wrong C uplo\n");
    if(d_C->nrows != d_A_left->nrows) eslog::error("incompatible matrices\n");

    n_size = n_end - n_start;
    k_start = h_A_pivots->vals[n_start];
    k_size = d_A_left->ncols - k_start;

    d_C_lower = d_C;
    if(d_C->prop.uplo == 'U') {
        d_C_reordered = d_C->get_transposed_reordered_view();
        d_C_lower = &d_C_reordered;
    }

    d_sub_C_herk.set_view(n_size, n_size, d_C_lower->ld, d_C_lower->order, nullptr);
    d_sub_C_herk.prop.uplo = d_C_lower->prop.uplo;
    d_sub_C_gemm.set_view(n_size, n_start, d_C_lower->ld, d_C_lower->order, nullptr);
    d_sub_A_left.set_view(n_size, k_size, d_A_left->ld, d_A_left->order, nullptr);
    d_sub_A_top.set_view(k_size, n_start, d_A_top->ld, d_A_top->order, nullptr);

    op_sub_C_herk.set_matrix_src(d_C_lower);
    op_sub_C_herk.set_matrix_dst(&d_sub_C_herk);
    op_sub_C_herk.set_bounds(n_start, n_end, n_start, n_end);

    op_sub_C_gemm.set_matrix_src(d_C_lower);
    op_sub_C_gemm.set_matrix_dst(&d_sub_C_gemm);
    op_sub_C_gemm.set_bounds(n_start, n_end, 0, n_start);

    op_sub_A_left.set_matrix_src(d_A_left);
    op_sub_A_left.set_matrix_dst(&d_sub_A_left);
    op_sub_A_left.set_bounds(n_start, n_end, k_start, d_A_left->ncols);

    op_sub_A_top.set_matrix_src(d_A_top);
    op_sub_A_top.set_matrix_dst(&d_sub_A_top);
    op_sub_A_top.set_bounds(k_start, d_A_top->nrows, 0, n_start);

    op_herk = herk_ddnx_ddny<T>::make();
    op_herk->set_handles(q, handle_dnblas);
    op_herk->set_matrix_A(&d_sub_A_left);
    op_herk->set_matrix_C(&d_sub_C_herk);
    op_herk->set_coefficients(alpha, beta);
    op_herk->set_mode(math::blas::herk_mode::AAh);
    op_herk->setup();
    wss_tmp_perform = std::max(wss_tmp_perform, op_herk->get_wss_tmp_perform());

    if(n_start > 0) {
        op_gemm = gemm_ddnx_ddny_ddnz<T>::make();
        op_gemm->set_handles(q, handle_dnblas);
        op_gemm->set_matrix_A(&d_sub_A_left);
        op_gemm->set_matrix_B(&d_sub_A_top);
        op_gemm->set_matrix_C(&d_sub_C_gemm);
        op_gemm->set_coefficients(alpha, beta);
        op_gemm->setup();
        wss_tmp_perform = std::max(wss_tmp_perform, op_gemm->get_wss_tmp_perform());
    }

    stacktimer::info("wss_tmp_perform    %zu", wss_tmp_perform);

    stacktimer::pop();

    called_setup = true;
}



template<typename T, typename I>
size_t gpu_herk_tri_chunk_stairs<T,I>::get_wss_tmp_perform()
{
    if(!called_setup) eslog::error("setup was not called\n");

    return wss_tmp_perform;
}



template<typename T, typename I>
void gpu_herk_tri_chunk_stairs<T,I>::perform_submit(void * ws_tmp)
{
    stacktimer::push("gpu_herk_tri_chunk_stairs::perform_submit");

    if(!called_setup) eslog::error("setup was not called\n");
    if(ws_tmp == nullptr && wss_tmp_perform > 0) eslog::error("temporary workspace is null\n");

    if(d_C->prop.uplo == 'U') {
        d_C_reordered = d_C->get_transposed_reordered_view();
    }

    op_sub_C_herk.perform();
    op_sub_C_gemm.perform();
    op_sub_A_left.perform();
    op_sub_A_top.perform();

    op_herk->perform_submit(ws_tmp);

    if(n_start > 0) {
        op_gemm->perform_submit(ws_tmp);
    }

    stacktimer::pop();
}



#define INSTANTIATE_T_I(T,I) \
template class gpu_herk_tri_chunk_stairs<T,I>;

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
