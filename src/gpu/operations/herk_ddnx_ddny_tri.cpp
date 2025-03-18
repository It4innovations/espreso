
#include "gpu/operations/herk_ddnx_ddny_tri.h"



namespace espreso {
namespace gpu {
namespace operations {



template<typename T, typename I>
void herk_ddnx_ddny_tri<T,I>::set_handles(gpu::mgm::queue q_, gpu::dnblas::handle dnblas_handle_)
{
    if(called_set_handles) eslog::error("handles are already set\n");

    q = q_;
    handle_dnblas = dnblas_handle_;

    called_set_handles = true;
}



template<typename T, typename I>
void herk_ddnx_ddny_tri<T,I>::set_matrix_A(MatrixDenseView_new<T> * A_)
{
    if(A != nullptr) eslog::error("matrix A is already set\n");
    if(A_ == nullptr) eslog::error("A cannot be nullptr\n");

    A = A_;
}



template<typename T, typename I>
void herk_ddnx_ddny_tri<T,I>::set_matrix_C(MatrixDenseView_new<T> * C_)
{
    if(C != nullptr) eslog::error("matrix C is already set\n");
    if(C_ == nullptr) eslog::error("C cannot be nullptr\n");

    C = C_;
}



template<typename T, typename I>
void herk_ddnx_ddny_tri<T,I>::set_A_pattern(MatrixCsxView<T,I> * A_pattern_host_)
{
    if(A_pattern_host != nullptr) eslog::error("A pattern is already set\n");
    if(A_pattern_host_ == nullptr) eslog::error("A pattern cannot be nullptr\n");

    A_pattern_host = A_pattern_host_;
}



template<typename T, typename I>
void herk_ddnx_ddny_tri<T,I>::set_coefficients(Treal alpha_, Treal beta_)
{
    alpha = alpha_;
    beta = beta_;
}



template<typename T, typename I>
void herk_ddnx_ddny_tri<T,I>::set_mode(math::herk_mode mode_)
{
    if(called_set_mode) eslog::error("mode is already set\n");

    mode = mode_;

    called_set_mode = true;
}



template<typename T, typename I>
void herk_ddnx_ddny_tri<T,I>::setup()
{
    if(!called_set_handles) eslog::error("handles are not set\n");
    if(A == nullptr) eslog::error("matrix A is not set\n");
    if(C == nullptr) eslog::error("matrix C is not set\n");
    if(A_pattern == nullptr) eslog::error("A pattern is not set\n");
    if(called_setup) eslog::error("setup was already called\n");
    if(h_L->nrows != h_L->ncols) eslog::error("matrix L is not square\n");
    if(h_L->nrows != d_X->nrows) eslog::error("incompatible matrices\n");
    if(h_L->prop.uplo != 'L') eslog::error("L has to have uplo=L\n");

    ator_ws_tmp_linear = std::make_unique<AllocatorArena_new>(false, true, gpu::mgm::get_natural_pitch_align());
    ator_ws_tmp_overlap = std::make_unique<AllocatorSinglePointer_new>(false, true, gpu::mgm::get_natural_pitch_align());

    A_reordered = A->get_transposed_reordered_view();

    if(mode == math::herk_mode::AAh) {
        A_left = A;
        A_right = A_reordered;
    }
    if(mode == math::herk_mode::AhA) {
        A_left = A_reordered;
        A_right = A;
    }

    n = A_left->nrows;
    k = A_left->ncols;
    
    A_pivots.set(A_pattern_host->ncols, AllocatorCPU_new::get_singleton());
    A_pivots.alloc();
    pivots_trails_csx<T,I>::do_all(A_top, &A_pivots, 'C', 'P', 'B');

    A_pivots.set(A_pattern_host->nrows, AllocatorCPU_new::get_singleton());
    A_pivots.alloc();
    pivots_trails_csx<T,I>::do_all(A_top, &A_pivots, 'R', 'T', 'F');

    for(size_t i = 1; i < A_pivots.size; i++) {
        if(A_pivots.vals[i-1] > A_pivots.vals[i]) {
            eslog::error("X does not have lower triangular structure\n");
        }
    }
    for(size_t i = 1; i < A_pivots.size; i++) {
        if(A_pivots.vals[i-1] > A_pivots.vals[i]) {
            eslog::error("X does not have lower triangular structure\n");
        }
    }

    tri_partition_herk partitioner;
    char partition_direction = '_';
    if(cfg.strategy == 'T') partition_direction = 'N';
    if(cfg.strategy == 'Q') partition_direction = 'K';
    partitioner.set_config(cfg.partition_algorithm, partition_direction, cfg.partition_parameter, cfg.strategy);
    partitioner.set_system(n, k);
    partitioner.set_output_partition(&partition);
    partitioner.setup();
    num_chunks = partitioner.get_num_chunks();
    partition.set(num_chunks + 1, AllocatorCPU_new::get_singleton());
    partition.alloc();
    partitioner.perform();

    if(cfg.strategy == 'T') {
        op_chunks_stairs.resize(num_chunks);
        for(size_t ch = 0; ch < num_chunks; ch++) {
            I n_start = partition.vals[ch];
            I n_end = partition.vals[ch+1];
            gpu_herk_tri_chunk_stairs<T,I> & op_chunk = op_chunks_stairs[ch];
            op_chunk.set_range(n_start, n_end);
            op_chunk.set_handles(q, handle_spblas);
            op_chunk.set_matrix_d_A_left(d_A_left);
            op_chunk.set_matrix_d_A_top(d_A_top);
            op_chunk.set_matrix_d_C(d_C);
            op_chunk.set_h_A_pivots(h_A_pivots);
            op_chunk.set_coefficients(alpha, beta);
            op_chunk.setup();
            wss_tmp_perform = std::max(wss_tmp_perform, op_chunk.get_wss_tmp_perform());
        }
    }
    if(cfg.strategy == 'Q') {
        d_A_top_dummy = d_A_top->get_submatrix_view(0, 0, 0, d_A_top->ncols);

        op_scale_C.set_handles(q, handle_dnblas);
        op_scale_C.set_matrix_A(&d_A_top_dummy);
        op_scale_C.set_matrix_C(d_C);
        op_scale_C.set_coefficients(Treal{0}, beta);
        op_scale_C.set_mode(math::herk_mode::AhA);
        op_scale_C.setup();
        wss_tmp_perform = std::max(wss_tmp_perform, op_scale_C.get_wss_tmp_perform());

        op_chunks_squares.resize(num_chunks);
        for(size_t ch = 0; ch < num_chunks; ch++) {
            I k_start = partition.vals[ch];
            I k_end = partition.vals[ch+1];
            gpu_herk_tri_chunk_squares<T,I> & op_chunk = op_chunks_squares[ch];
            op_chunk.set_range(k_start, k_end);
            op_chunk.set_handles(q, handle_spblas);
            op_chunk.set_matrix_d_A_left(d_A_left);
            op_chunk.set_matrix_d_A_top(d_A_top);
            op_chunk.set_matrix_d_C(d_C);
            op_chunk.set_h_A_trails(h_A_trails);
            op_chunk.set_coefficients(alpha);
            op_chunk.setup();
            wss_tmp_perform = std::max(wss_tmp_perform, op_chunk.get_wss_tmp_perform());
        }
    }

    wss_tmp_perform_linear = ((wss_tmp_perform_linear - 1) / ator_ws_tmp_linear->get_align() + 1) * ator_ws_tmp_linear->get_align();
    wss_tmp_perform = wss_tmp_perform_linear + wss_tmp_perform_overlap;

    called_setup = true;
}



template<typename T, typename I>
size_t herk_ddnx_ddny_tri<T,I>::get_wss_tmp_perform()
{
    if(!called_setup) eslog::error("setup was not called\n");

    return wss_tmp_perform;
}



template<typename T, typename I>
void herk_ddnx_ddny_tri<T,I>::perform_submit(void * ws_tmp)
{
    if(!called_preprocess) eslog::error("preprocess has not been called\n");
    if(ws_tmp == nullptr && wss_tmp_perform > 0) eslog::error("temporary workspace is null\n");

    ator_ws_tmp_linear.set(ws_tmp, wss_tmp_perform_linear);
    ator_ws_tmp_overlap.set(ws_tmp + wss_tmp_perform_linear, wss_tmp_perform_overlap);

    A_reordered = A->get_transposed_reordered_view();

    if(cfg.strategy == 'T') {
        for(size_t ch = 0; ch < num_chunks; ch++) {
            op_chunks_stairs[ch].perform_submit(ws_tmp);
        }
    }
    if(cfg.strategy == 'Q') {
        d_A_top_dummy = d_A_top->get_submatrix_view(0, 0, 0, d_A_top->ncols);
        op_scale_C.perform_submit(ws_tmp);
        for(size_t ch = 0; ch < num_chunks; ch++) {
            op_chunks_squares[ch].perform_submit(ws_tmp);
        }
    }

    ator_ws_tmp_linear.unset();
    ator_ws_tmp_overlap.unset();
}



#define INSTANTIATE_T_I(T,I) \
template class herk_ddnx_ddny_tri<T,I>;

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

#endif /* SRC_GPU_OPERATIONS_HERK_DDNX_DDNY_H */
