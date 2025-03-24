
#include "gpu/operations/auxiliary/gpu_herk_tri_chunk_squares.h"

#include "basis/utilities/stacktimer.h"



namespace espreso {
namespace gpu {
namespace operations {



template<typename T, typename I>
void gpu_herk_tri_chunk_squares<T,I>::set_range(size_t k_start_, size_t k_end_)
{
    if(called_set_range) eslog::error("range is already set\n");

    k_start = k_start_;
    k_end = k_end_;

    called_set_range = true;
}



template<typename T, typename I>
void gpu_herk_tri_chunk_squares<T,I>::set_handles(gpu::mgm::queue q_, gpu::dnblas::handle handle_dnblas_)
{
    if(called_set_handles) eslog::error("handles are already set\n");

    q = q_;
    handle_dnblas = handle_dnblas_;

    called_set_handles = true;
}



template<typename T, typename I>
void gpu_herk_tri_chunk_squares<T,I>::set_matrix_d_A_left(MatrixDenseView_new<T> * d_A_left_)
{
    if(d_A_left != nullptr) eslog::error("d_A_left is already set\n");
    if(d_A_left_ == nullptr) eslog::error("d_A_left cannot be nullptr\n");

    d_A_left = d_A_left_;
}



template<typename T, typename I>
void gpu_herk_tri_chunk_squares<T,I>::set_matrix_d_A_top(MatrixDenseView_new<T> * d_A_top_)
{
    if(d_A_top != nullptr) eslog::error("d_A_top is already set\n");
    if(d_A_top_ == nullptr) eslog::error("d_A_top cannot be nullptr\n");

    d_A_top = d_A_top_;
}



template<typename T, typename I>
void gpu_herk_tri_chunk_squares<T,I>::set_matrix_d_C(MatrixDenseView_new<T> * d_C_)
{
    if(d_C != nullptr) eslog::error("d_C is already set\n");
    if(d_C_ == nullptr) eslog::error("d_C cannot be nullptr\n");

    d_C = d_C_;
}



template<typename T, typename I>
void gpu_herk_tri_chunk_squares<T,I>::set_h_A_trails(VectorDenseView_new<I> * h_A_trails_)
{
    if(h_A_trails != nullptr) eslog::error("h_A_trails is already set\n");
    if(h_A_trails_ == nullptr) eslog::error("h_A_trails cannot be nullptr\n");

    h_A_trails = h_A_trails_;
}



template<typename T, typename I>
void gpu_herk_tri_chunk_squares<T,I>::set_coefficients(Treal alpha_)
{
    alpha = alpha_;
}



template<typename T, typename I>
void gpu_herk_tri_chunk_squares<T,I>::setup()
{
    stacktimer::push("gpu_herk_tri_chunk_squares::setup");

    if(utils::is_complex<T>()) eslog::error("complex numbers not supported\n");
    if(!called_set_range) eslog::error("range is not set\n");
    if(!called_set_handles) eslog::error("handles are not set\n");
    if(called_setup) eslog::error("setup was already called\n");
    if(d_A_left == nullptr) eslog::error("d_A_left is not set\n");
    if(d_A_top == nullptr) eslog::error("d_A_top is not set\n");
    if(d_C == nullptr) eslog::error("d_C is not set\n");
    if(h_A_trails == nullptr) eslog::error("h_A_trails is not set\n");
    if(d_A_left->nrows != d_A_top->ncols || d_A_left->ncols != d_A_top->nrows) eslog::error("wrong matrices A\n");
    if(d_A_left->order == d_A_top->order) eslog::error("wrong matrices A\n");
    if(d_A_left->vals != d_A_top->vals) eslog::error("wrong matrices A\n");
    if(d_C->nrows != d_C->ncols) eslog::error("d_C must be square\n");
    if(d_C->prop.uplo != 'U' && d_C->prop.uplo != 'L') eslog::error("wrong C uplo\n");
    if(d_C->nrows != d_A_left->nrows) eslog::error("incompatible matrices\n");

    k_size = k_end - k_start;
    n_end = h_A_trails->vals[k_end - 1] + 1;
    n_size = n_end;

    d_sub_C.set_view(n_size, n_size, d_C->ld, d_C->order, nullptr);
    d_sub_C.prop.uplo = d_C->prop.uplo;
    d_sub_A_top.set_view(k_size, n_size, d_A_top->ld, d_A_top->order, nullptr);

    op_sub_C.set_matrix_src(d_C);
    op_sub_C.set_matrix_dst(&d_sub_C);
    op_sub_C.set_bounds(0, n_end, 0, n_end);

    op_sub_A_top.set_matrix_src(d_A_top);
    op_sub_A_top.set_matrix_dst(&d_sub_A_top);
    op_sub_A_top.set_bounds(k_start, k_end, 0, n_end);

    op_herk = herk_ddnx_ddny<T>::make();
    op_herk->set_handles(q, handle_dnblas);
    op_herk->set_matrix_A(&d_sub_A_top);
    op_herk->set_matrix_C(&d_sub_C);
    op_herk->set_coefficients(alpha, Treal{1});
    op_herk->set_mode(math::blas::herk_mode::AhA);
    op_herk->setup();
    wss_tmp_perform = std::max(wss_tmp_perform, op_herk->get_wss_tmp_perform());

    stacktimer::pop();

    called_setup = true;
}



template<typename T, typename I>
size_t gpu_herk_tri_chunk_squares<T,I>::get_wss_tmp_perform()
{
    if(!called_setup) eslog::error("setup was not called\n");

    return wss_tmp_perform;
}



template<typename T, typename I>
void gpu_herk_tri_chunk_squares<T,I>::perform_submit(void * ws_tmp)
{
    stacktimer::push("gpu_herk_tri_chunk_squares::perform_submit");

    if(!called_setup) eslog::error("setup was not called\n");
    if(ws_tmp == nullptr && wss_tmp_perform > 0) eslog::error("temporary workspace is null\n");

    op_sub_C.perform();
    op_sub_A_top.perform();

    op_herk->perform_submit(ws_tmp);

    stacktimer::pop();
}



#define INSTANTIATE_T_I(T,I) \
template class gpu_herk_tri_chunk_squares<T,I>;

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
