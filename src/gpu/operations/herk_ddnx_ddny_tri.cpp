
#include "gpu/operations/herk_ddnx_ddny_tri.h"

#include "basis/utilities/stacktimer.h"
#include "math/operations/pivots_trails_csx.h"
#include "math/operations/auxiliary/tri_partition_herk.h"



namespace espreso {
namespace gpu {
namespace operations {



template<typename T, typename I>
void herk_ddnx_ddny_tri<T,I>::set_config(config cfg_)
{
    if(called_set_config) eslog::error("config is already set\n");

    cfg = cfg_;

    called_set_config = true;
}



template<typename T, typename I>
void herk_ddnx_ddny_tri<T,I>::set_handles(gpu::mgm::queue q_, gpu::dnblas::handle dnblas_handle_)
{
    if(called_set_handles) eslog::error("handles are already set\n");

    q = q_;
    handle_dnblas = dnblas_handle_;

    called_set_handles = true;
}



template<typename T, typename I>
void herk_ddnx_ddny_tri<T,I>::set_matrix_d_A(MatrixDenseView_new<T> * d_A_)
{
    if(d_A != nullptr) eslog::error("matrix d_A is already set\n");
    if(d_A_ == nullptr) eslog::error("d_A cannot be nullptr\n");

    d_A = d_A_;
}



template<typename T, typename I>
void herk_ddnx_ddny_tri<T,I>::set_matrix_d_C(MatrixDenseView_new<T> * d_C_)
{
    if(d_C != nullptr) eslog::error("matrix d_C is already set\n");
    if(d_C_ == nullptr) eslog::error("d_C cannot be nullptr\n");

    d_C = d_C_;
}



template<typename T, typename I>
void herk_ddnx_ddny_tri<T,I>::set_h_A_pattern(MatrixCsxView_new<T,I> * h_A_pattern_)
{
    if(h_A_pattern != nullptr) eslog::error("A pattern is already set\n");
    if(h_A_pattern_ == nullptr) eslog::error("A pattern cannot be nullptr\n");

    h_A_pattern = h_A_pattern_;
}



template<typename T, typename I>
void herk_ddnx_ddny_tri<T,I>::set_coefficients(Treal alpha_, Treal beta_)
{
    alpha = alpha_;
    beta = beta_;
}



template<typename T, typename I>
void herk_ddnx_ddny_tri<T,I>::set_mode(math::blas::herk_mode mode_)
{
    if(called_set_mode) eslog::error("mode is already set\n");

    mode = mode_;

    called_set_mode = true;
}



template<typename T, typename I>
void herk_ddnx_ddny_tri<T,I>::setup()
{
    stacktimer::push("herk_ddnx_ddny_tri::setup");

    if(!called_set_config) eslog::error("config is not set\n");
    if(!called_set_handles) eslog::error("handles are not set\n");
    if(d_A == nullptr) eslog::error("matrix A is not set\n");
    if(d_C == nullptr) eslog::error("matrix C is not set\n");
    if(h_A_pattern == nullptr) eslog::error("A pattern is not set\n");
    if(!d_A->ator->is_data_accessible_gpu()) eslog::error("matrix d_A must be gpu-accessible\n");
    if(!d_C->ator->is_data_accessible_gpu()) eslog::error("matrix d_C must be gpu-accessible\n");
    if(!h_A_pattern->ator->is_data_accessible_cpu()) eslog::error("h_A_pattern must be cpu-accessible\n");
    if(called_setup) eslog::error("setup was already called\n");
    if(d_C->nrows != d_C->ncols) eslog::error("matrix C is not square\n");
    if(mode == math::blas::herk_mode::AhA && d_A->ncols != d_C->ncols) eslog::error("incompatible matrices\n");
    if(mode == math::blas::herk_mode::AAh && d_A->nrows != d_C->nrows) eslog::error("incompatible matrices\n");
    if(d_C->prop.uplo != 'L' && d_C->prop.uplo != 'U') eslog::error("invalid matrix C uplo\n");

    ator_ws_tmp_linear = std::make_unique<AllocatorArena_new>(AllocatorGPU_new::get_singleton());
    ator_ws_tmp_overlap = std::make_unique<AllocatorSinglePointer_new>(AllocatorGPU_new::get_singleton());

    d_A_reordered = d_A->get_transposed_reordered_view();

    if(mode == math::blas::herk_mode::AAh) {
        d_A_left = d_A;
        d_A_top = &d_A_reordered;
    }
    if(mode == math::blas::herk_mode::AhA) {
        d_A_left = &d_A_reordered;
        d_A_top = d_A;
    }

    n = d_A_left->nrows;
    k = d_A_left->ncols;
    
    h_A_pivots.set(h_A_pattern->ncols, AllocatorCPU_new::get_singleton());
    h_A_pivots.alloc();
    if(mode == math::blas::herk_mode::AhA) math::operations::pivots_trails_csx<T,I>::do_all(h_A_pattern, &h_A_pivots, 'C', 'P', 'B');
    if(mode == math::blas::herk_mode::AAh) math::operations::pivots_trails_csx<T,I>::do_all(h_A_pattern, &h_A_pivots, 'R', 'P', 'B');

    h_A_trails.set(h_A_pattern->nrows, AllocatorCPU_new::get_singleton());
    h_A_trails.alloc();
    if(mode == math::blas::herk_mode::AhA) math::operations::pivots_trails_csx<T,I>::do_all(h_A_pattern, &h_A_trails, 'R', 'T', 'F');
    if(mode == math::blas::herk_mode::AAh) math::operations::pivots_trails_csx<T,I>::do_all(h_A_pattern, &h_A_trails, 'C', 'T', 'F');

    for(size_t i = 1; i < h_A_pivots.size; i++) {
        if(h_A_pivots.vals[i-1] > h_A_pivots.vals[i]) {
            eslog::error("A does not have correct triangular structure\n");
        }
    }
    for(size_t i = 1; i < h_A_trails.size; i++) {
        if(h_A_trails.vals[i-1] > h_A_trails.vals[i]) {
            eslog::error("A does not have correct triangular structure\n");
        }
    }

    math::operations::tri_partition_herk partitioner;
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
            stacktimer::info("herk_ddnx_ddny_tri::setup chunk %zu", ch);
    
            I n_start = partition.vals[ch];
            I n_end = partition.vals[ch+1];
            gpu_herk_tri_chunk_stairs<T,I> & op_chunk = op_chunks_stairs[ch];
            op_chunk.set_range(n_start, n_end);
            op_chunk.set_handles(q, handle_dnblas);
            op_chunk.set_matrix_d_A_left(d_A_left);
            op_chunk.set_matrix_d_A_top(d_A_top);
            op_chunk.set_matrix_d_C(d_C);
            op_chunk.set_h_A_pivots(&h_A_pivots);
            op_chunk.set_coefficients(alpha, beta);
            op_chunk.setup();
            wss_tmp_perform = std::max(wss_tmp_perform, op_chunk.get_wss_tmp_perform());
        }
    }
    if(cfg.strategy == 'Q') {
        d_A_top_dummy = d_A_top->get_submatrix_view(0, 0, 0, d_A_top->ncols);

        op_scale_C = herk_ddnx_ddny<T>::make();
        op_scale_C->set_handles(q, handle_dnblas);
        op_scale_C->set_matrix_A(&d_A_top_dummy);
        op_scale_C->set_matrix_C(d_C);
        op_scale_C->set_coefficients(Treal{0}, beta);
        op_scale_C->set_mode(math::blas::herk_mode::AhA);
        op_scale_C->setup();
        wss_tmp_perform = std::max(wss_tmp_perform, op_scale_C->get_wss_tmp_perform());

        op_chunks_squares.resize(num_chunks);
        for(size_t ch = 0; ch < num_chunks; ch++) {
            stacktimer::info("herk_ddnx_ddny_tri::setup chunk %zu", ch);
    
            I k_start = partition.vals[ch];
            I k_end = partition.vals[ch+1];
            gpu_herk_tri_chunk_squares<T,I> & op_chunk = op_chunks_squares[ch];
            op_chunk.set_range(k_start, k_end);
            op_chunk.set_handles(q, handle_dnblas);
            op_chunk.set_matrix_d_A_left(d_A_left);
            op_chunk.set_matrix_d_A_top(d_A_top);
            op_chunk.set_matrix_d_C(d_C);
            op_chunk.set_h_A_trails(&h_A_trails);
            op_chunk.set_coefficients(alpha);
            op_chunk.setup();
            wss_tmp_perform = std::max(wss_tmp_perform, op_chunk.get_wss_tmp_perform());
        }
    }

    wss_tmp_perform_linear = utils::round_up(wss_tmp_perform_linear, ator_ws_tmp_linear->get_align());
    wss_tmp_perform = wss_tmp_perform_linear + wss_tmp_perform_overlap;

    stacktimer::info("wss_tmp_perform    %zu", wss_tmp_perform);

    stacktimer::pop();

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
    stacktimer::push("herk_ddnx_ddny_tri::perform_submit");

    if(!called_setup) eslog::error("setup has not been called\n");
    if(ws_tmp == nullptr && wss_tmp_perform > 0) eslog::error("temporary workspace is null\n");

    ator_ws_tmp_linear->set(ws_tmp, wss_tmp_perform_linear);
    ator_ws_tmp_overlap->set((char*)ws_tmp + wss_tmp_perform_linear, wss_tmp_perform_overlap);

    d_A_reordered = d_A->get_transposed_reordered_view();

    if(cfg.strategy == 'T') {
        for(size_t ch = 0; ch < num_chunks; ch++) {
            stacktimer::info("herk_ddnx_ddny_tri::perform_submit chunk %zu", ch);
    
            op_chunks_stairs[ch].perform_submit(ws_tmp);
        }
    }
    if(cfg.strategy == 'Q') {
        d_A_top_dummy = d_A_top->get_submatrix_view(0, 0, 0, d_A_top->ncols);
        op_scale_C->perform_submit(ws_tmp);
        for(size_t ch = 0; ch < num_chunks; ch++) {
            stacktimer::info("herk_ddnx_ddny_tri::perform_submit chunk %zu", ch);

            op_chunks_squares[ch].perform_submit(ws_tmp);
        }
    }

    ator_ws_tmp_linear->unset();
    ator_ws_tmp_overlap->unset();

    stacktimer::pop();
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
        INSTANTIATE_T(std::complex<double>)

            INSTANTIATE

        #undef INSTANTIATE
    #undef INSTANTIATE_T
#undef INSTANTIATE_T_I
    


}
}
}
