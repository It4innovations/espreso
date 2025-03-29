
#include "gpu/operations/auxiliary/trsm_hcsx_ddny_tri_splitfactor.h"

#include "basis/utilities/stacktimer.h"
#include "math/operations/pivots_trails_csx.h"
#include "math/operations/auxiliary/tri_partition_trsm.h"
#include "math/operations/submatrix_csx_csy.h"
#include "math/operations/pruning_subset_csx.h"



namespace espreso {
namespace gpu {
namespace operations {



template<typename T, typename I>
void trsm_hcsx_ddny_tri_splitfactor<T,I>::set_config(config cfg_)
{
    cfg = cfg_;

    called_set_config = true;
}



template<typename T, typename I>
void trsm_hcsx_ddny_tri_splitfactor<T,I>::set_handles(gpu::mgm::queue q_, gpu::spblas::handle spblas_handle_, gpu::dnblas::handle dnblas_handle_)
{
    if(called_set_handles) eslog::error("handles are already set\n");

    q = q_;
    handle_spblas = spblas_handle_;
    handle_dnblas = dnblas_handle_;

    called_set_handles = true;
}



template<typename T, typename I>
void trsm_hcsx_ddny_tri_splitfactor<T,I>::set_matrix_h_L(MatrixCsxView_new<T,I> * h_L_)
{
    if(h_L != nullptr) eslog::error("matrix h_L is already set\n");
    if(h_L_ == nullptr) eslog::error("h_L cannot be nullptr\n");

    h_L = h_L_;
}



template<typename T, typename I>
void trsm_hcsx_ddny_tri_splitfactor<T,I>::set_matrix_d_X(MatrixDenseView_new<T> * d_X_)
{
    if(d_X != nullptr) eslog::error("matrix d_X is already set\n");
    if(d_X_ == nullptr) eslog::error("d_X cannot be nullptr\n");

    d_X = d_X_;
}



template<typename T, typename I>
void trsm_hcsx_ddny_tri_splitfactor<T,I>::set_h_X_pattern(MatrixCsxView_new<T,I> * h_X_pattern_)
{
    if(h_X_pattern != nullptr) eslog::error("matrix h_X_pattern is already set\n");
    if(h_X_pattern_ == nullptr) eslog::error("h_X_pattern cannot be nullptr\n");

    h_X_pattern = h_X_pattern_;
}



template<typename T, typename I>
void trsm_hcsx_ddny_tri_splitfactor<T,I>::setup()
{
    stacktimer::push("trsm_hcsx_ddny_tri_splitfactor::setup");

    if(!called_set_handles) eslog::error("handles are not set\n");
    if(h_L == nullptr) eslog::error("matrix L is not set\n");
    if(d_X == nullptr) eslog::error("matrix X is not set\n");
    if(h_X_pattern == nullptr) eslog::error("X pattern is not set\n");
    if(called_setup) eslog::error("setup was already called\n");
    if(h_L->nrows != h_L->ncols) eslog::error("matrix L is not square\n");
    if(h_L->nrows != d_X->nrows) eslog::error("incompatible matrices\n");
    if(h_L->prop.uplo != 'L') eslog::error("L has to have uplo=L\n");

    ator_ws_persistent = std::make_unique<AllocatorArena_new>(false, true, gpu::mgm::get_natural_pitch_align());
    ator_ws_tmp_linear = std::make_unique<AllocatorArena_new>(false, true, gpu::mgm::get_natural_pitch_align());
    ator_ws_tmp_overlap = std::make_unique<AllocatorSinglePointer_new>(false, true, gpu::mgm::get_natural_pitch_align());

    h_X_colpivots.set(h_X_pattern->ncols, AllocatorCPU_new::get_singleton());
    h_X_colpivots.alloc();
    math::operations::pivots_trails_csx<T,I>::do_all(h_X_pattern, &h_X_colpivots, 'C', 'P', 'B');

    h_X_rowtrails.set(h_X_pattern->nrows, AllocatorCPU_new::get_singleton());
    h_X_rowtrails.alloc();
    math::operations::pivots_trails_csx<T,I>::do_all(h_X_pattern, &h_X_rowtrails, 'R', 'T', 'F');

    for(size_t i = 1; i < h_X_colpivots.size; i++) {
        if(h_X_colpivots.vals[i-1] > h_X_colpivots.vals[i]) {
            eslog::error("X does not have lower triangular structure\n");
        }
    }
    for(size_t i = 1; i < h_X_rowtrails.size; i++) {
        if(h_X_rowtrails.vals[i-1] > h_X_rowtrails.vals[i]) {
            eslog::error("X does not have lower triangular structure\n");
        }
    }

    math::operations::tri_partition_trsm partitioner;
    partitioner.set_config(cfg.partition.algorithm, 'V', cfg.partition.parameter);
    partitioner.set_system(d_X->nrows, d_X->ncols);
    partitioner.set_output_partition(&partition);
    partitioner.setup();
    num_chunks = partitioner.get_num_chunks();
    partition.set(num_chunks + 1, AllocatorCPU_new::get_singleton());
    partition.alloc();
    partitioner.perform();

    ops_chunks.resize(num_chunks);

    for(size_t ch = 0; ch < num_chunks; ch++) {
        stacktimer::info("trsm_hcsx_ddny_tri_splitfactor::setup chunk %zu", ch);

        size_t k_start = partition.vals[ch];
        size_t k_end = partition.vals[ch+1];

        typename gpu_trsm_trirhs_chunk_splitfactor<T,I>::config op_config;
        op_config.trsm_factor_spdn = cfg.trsm_factor_spdn;
        op_config.trsm_factor_order = cfg.trsm_factor_order;
        op_config.gemm_factor_prune = cfg.gemm_factor_prune;
        
        if(cfg.gemm_spdn_criteria == 'S') {
            op_config.gemm_factor_spdn = 'S';
        }
        if(cfg.gemm_spdn_criteria == 'D') {
            op_config.gemm_factor_spdn = 'D';
        }
        if(cfg.gemm_spdn_criteria == 'C') {
            double fraction_treshold = cfg.gemm_spdn_param;
            double curr_fraction = ((double)ch + 0.5) / num_chunks;
            op_config.gemm_factor_spdn = ((curr_fraction < fraction_treshold) ? 'S' : 'D');
        }
        if(cfg.gemm_spdn_criteria == 'Z') {
            double fraction_treshold = cfg.gemm_spdn_param;
            double curr_fraction = (double)(k_start + k_end) / 2.0 / d_X->nrows;
            op_config.gemm_factor_spdn = ((curr_fraction < fraction_treshold) ? 'S' : 'D');
        }
        if(cfg.gemm_spdn_criteria == 'T') {
            math::operations::submatrix_csx_csy<T,I> op_sub_L_bot;
            op_sub_L_bot.set_matrix_src(h_L);
            op_sub_L_bot.set_bounds(k_end, h_L->nrows, k_start, k_end);
            op_sub_L_bot.setup();
            size_t nnz = op_sub_L_bot.get_output_matrix_nnz();
            if(cfg.gemm_factor_prune == 'N') {
                size_t nvals = (h_L->nrows - k_end) * (k_end - k_start);
                double fraction_treshold = cfg.gemm_spdn_param;
                double curr_fraction = (double)nnz / nvals;
                op_config.gemm_factor_spdn = ((curr_fraction < fraction_treshold) ? 'S' : 'D');
            }
            if(cfg.gemm_factor_prune != 'N') {
                MatrixCsxData_new<T,I> sub_L_bot_test;
                sub_L_bot_test.set(h_L->nrows - k_end, k_end - k_start, nnz, 'R', AllocatorCPU_new::get_singleton());
                sub_L_bot_test.alloc();
                op_sub_L_bot.set_matrix_dst(&sub_L_bot_test);
                op_sub_L_bot.perform();
                math::operations::pruning_subset_csx<T,I> op_pruning_subset;
                op_pruning_subset.set_matrix(&sub_L_bot_test);
                char pm = cfg.gemm_factor_prune;
                op_pruning_subset.set_pruning_mode(pm == 'R' || pm == 'A', pm == 'C' || pm == 'A');
                op_pruning_subset.setup();
                size_t nvals = op_pruning_subset.get_pruned_nrows() * op_pruning_subset.get_pruned_ncols();
                sub_L_bot_test.clear();
                double fraction_treshold = cfg.gemm_spdn_param;
                double curr_fraction = (double)nnz / nvals;
                op_config.gemm_factor_spdn = ((curr_fraction < fraction_treshold) ? 'S' : 'D');
            }
        }
        
        if(op_config.gemm_factor_spdn == 'S') {
            op_config.gemm_factor_order = cfg.gemm_factor_order_sp;
        }
        if(op_config.gemm_factor_spdn == 'D') {
            op_config.gemm_factor_order = cfg.gemm_factor_order_dn;
        }

        gpu_trsm_trirhs_chunk_splitfactor<T,I> & op_chunk = ops_chunks[ch];
        op_chunk.set_config(op_config);
        op_chunk.set_range(k_start, k_end);
        op_chunk.set_handles(q, handle_spblas, handle_dnblas);
        op_chunk.set_matrix_h_L(h_L);
        op_chunk.set_matrix_d_X(d_X);
        op_chunk.set_h_X_rowtrails(&h_X_rowtrails);
        op_chunk.setup();
        wss_internal += op_chunk.get_wss_internal();
        wss_persistent += utils::round_up(op_chunk.get_wss_persistent(), ator_ws_persistent->get_align());
        wss_tmp_preprocess_overlap = std::max(wss_tmp_preprocess_overlap, op_chunk.get_wss_tmp_preprocess());
        wss_tmp_perform_overlap = std::max(wss_tmp_perform_overlap, op_chunk.get_wss_tmp_perform());
    }

    wss_tmp_preprocess_linear = utils::round_up(wss_tmp_preprocess_linear, ator_ws_tmp_linear->get_align());
    wss_tmp_perform_linear = utils::round_up(wss_tmp_perform_linear, ator_ws_tmp_linear->get_align());

    wss_tmp_preprocess = wss_tmp_preprocess_linear + wss_tmp_preprocess_overlap;
    wss_tmp_perform = wss_tmp_perform_linear + wss_tmp_perform_overlap;

    stacktimer::info("wss_internal       %zu", wss_internal);
    stacktimer::info("wss_persistent     %zu", wss_persistent);
    stacktimer::info("wss_tmp_preprocess %zu", wss_tmp_preprocess);
    stacktimer::info("wss_tmp_perform    %zu", wss_tmp_perform);

    stacktimer::pop();

    called_setup = true;
}



template<typename T, typename I>
size_t trsm_hcsx_ddny_tri_splitfactor<T,I>::get_wss_internal()
{
    if(!called_setup) eslog::error("setup was not called\n");

    return wss_internal;
}



template<typename T, typename I>
size_t trsm_hcsx_ddny_tri_splitfactor<T,I>::get_wss_persistent()
{
    if(!called_setup) eslog::error("setup was not called\n");

    return wss_persistent;
}



template<typename T, typename I>
size_t trsm_hcsx_ddny_tri_splitfactor<T,I>::get_wss_tmp_preprocess()
{
    if(!called_setup) eslog::error("setup was not called\n");

    return wss_tmp_preprocess;
}



template<typename T, typename I>
size_t trsm_hcsx_ddny_tri_splitfactor<T,I>::get_wss_tmp_perform()
{
    if(!called_setup) eslog::error("setup was not called\n");

    return wss_tmp_perform;
}



template<typename T, typename I>
void trsm_hcsx_ddny_tri_splitfactor<T,I>::set_ws_persistent(void * ws_persistent_)
{
    if(ws_persistent_ == nullptr && wss_persistent > 0) eslog::error("persistent workspace is null\n");
    if(ws_persistent != nullptr) eslog::error("cannot re-set persistent workspace\n");

    ws_persistent = ws_persistent_;
}



template<typename T, typename I>
void trsm_hcsx_ddny_tri_splitfactor<T,I>::preprocess_submit(void * ws_tmp)
{
    stacktimer::push("trsm_hcsx_ddny_tri_splitfactor::preprocess_submit");

    if(!called_setup) eslog::error("setup has not been called\n");
    if(called_preprocess) eslog::error("preprocess has already been called\n");
    if(ws_tmp == nullptr && wss_tmp_preprocess > 0) eslog::error("temporary workspace is null\n");

    ator_ws_persistent->set(ws_persistent, wss_persistent);

    ator_ws_tmp_linear->set(ws_tmp, wss_tmp_preprocess_linear);
    ator_ws_tmp_overlap->set((char*)ws_tmp + wss_tmp_preprocess_linear, wss_tmp_preprocess_overlap);

    for(size_t ch = 0; ch < num_chunks; ch++) {
        stacktimer::info("trsm_hcsx_ddny_tri_splitfactor::preprocess_submit chunk %zu", ch);

        gpu_trsm_trirhs_chunk_splitfactor<T,I> & op_chunk = ops_chunks[ch];
        op_chunk.set_ws_persistent(ator_ws_persistent->alloc(op_chunk.get_wss_persistent()));
        op_chunk.preprocess_submit(ator_ws_tmp_overlap->alloc(op_chunk.get_wss_tmp_preprocess()));
    }
    
    ator_ws_tmp_linear->unset();
    ator_ws_tmp_overlap->unset();

    stacktimer::pop();

    called_preprocess = true;
}



template<typename T, typename I>
void trsm_hcsx_ddny_tri_splitfactor<T,I>::perform_submit(void * ws_tmp)
{
    stacktimer::push("trsm_hcsx_ddny_tri_splitfactor::perform_submit");

    if(!called_preprocess) eslog::error("preprocess has not been called\n");
    if(ws_tmp == nullptr && wss_tmp_perform > 0) eslog::error("temporary workspace is null\n");

    ator_ws_tmp_linear->set(ws_tmp, wss_tmp_perform_linear);
    ator_ws_tmp_overlap->set((char*)ws_tmp + wss_tmp_perform_linear, wss_tmp_perform_overlap);

    for(size_t ch = 0; ch < num_chunks; ch++) {
        stacktimer::info("trsm_hcsx_ddny_tri_splitfactor::perform_submit chunk %zu", ch);

        gpu_trsm_trirhs_chunk_splitfactor<T,I> & op_chunk = ops_chunks[ch];
        op_chunk.perform_submit(ator_ws_tmp_overlap->alloc(op_chunk.get_wss_tmp_perform()));
    }

    ator_ws_tmp_linear->unset();
    ator_ws_tmp_overlap->unset();

    stacktimer::pop();
}



#define INSTANTIATE_T_I(T,I) \
template class trsm_hcsx_ddny_tri_splitfactor<T,I>;

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
