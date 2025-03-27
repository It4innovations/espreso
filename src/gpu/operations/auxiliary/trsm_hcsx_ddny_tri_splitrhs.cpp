
#include "gpu/operations/auxiliary/trsm_hcsx_ddny_tri_splitrhs.h"

#include "basis/utilities/stacktimer.h"
#include "math/operations/pivots_trails_csx.h"
#include "math/operations/auxiliary/tri_partition_trsm.h"
#include "math/operations/submatrix_csx_csy.h"



namespace espreso {
namespace gpu {
namespace operations {



template<typename T, typename I>
void trsm_hcsx_ddny_tri_splitrhs<T,I>::set_config(config cfg_)
{
    cfg = cfg_;

    called_set_config = true;
}



template<typename T, typename I>
void trsm_hcsx_ddny_tri_splitrhs<T,I>::set_handles(gpu::mgm::queue q_, gpu::spblas::handle spblas_handle_, gpu::dnblas::handle dnblas_handle_)
{
    if(called_set_handles) eslog::error("handles are already set\n");

    q = q_;
    handle_spblas = spblas_handle_;
    handle_dnblas = dnblas_handle_;

    called_set_handles = true;
}



template<typename T, typename I>
void trsm_hcsx_ddny_tri_splitrhs<T,I>::set_matrix_h_L(MatrixCsxView_new<T,I> * h_L_)
{
    if(h_L != nullptr) eslog::error("matrix h_L is already set\n");
    if(h_L_ == nullptr) eslog::error("h_L cannot be nullptr\n");

    h_L = h_L_;
}



template<typename T, typename I>
void trsm_hcsx_ddny_tri_splitrhs<T,I>::set_matrix_d_X(MatrixDenseView_new<T> * d_X_)
{
    if(d_X != nullptr) eslog::error("matrix d_X is already set\n");
    if(d_X_ == nullptr) eslog::error("d_X cannot be nullptr\n");

    d_X = d_X_;
}



template<typename T, typename I>
void trsm_hcsx_ddny_tri_splitrhs<T,I>::set_h_X_pattern(MatrixCsxView_new<T,I> * h_X_pattern_)
{
    if(h_X_pattern != nullptr) eslog::error("matrix h_X_pattern is already set\n");
    if(h_X_pattern_ == nullptr) eslog::error("h_X_pattern cannot be nullptr\n");

    h_X_pattern = h_X_pattern_;
}



template<typename T, typename I>
void trsm_hcsx_ddny_tri_splitrhs<T,I>::setup()
{
    stacktimer::push("trsm_hcsx_ddny_tri_splitrhs::setup");

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
    partitioner.set_config(cfg.partition.algorithm, 'H', cfg.partition.parameter);
    partitioner.set_system(d_X->nrows, d_X->ncols);
    partitioner.set_output_partition(&partition);
    partitioner.setup();
    num_chunks = partitioner.get_num_chunks();
    partition.set(num_chunks + 1, AllocatorCPU_new::get_singleton());
    partition.alloc();
    partitioner.perform();

    h_L_to_use = h_L;

    if(cfg.factor_order_sp != h_L->order || cfg.factor_order_dn != h_L->order) {
        h_L_reordered.set(h_L->nrows, h_L->ncols, h_L->nnz, change_order(h_L->order), AllocatorCPU_new::get_singleton());
        h_L_reordered.prop.uplo = h_L->prop.uplo;
        h_L_reordered.prop.diag = h_L->prop.diag;
        h_L_reordered.alloc();
        op_h_L_reorder.set_matrix_src(h_L);
        op_h_L_reorder.set_matrix_dst(&h_L_reordered);
        op_h_L_reorder.perform_pattern();
    }
    if(cfg.factor_order_sp != h_L->order) {
        h_L_to_use = &h_L_reordered;
    }
    
    h_L_nnzinsubs_splitrhs.set(h_L->ncols, AllocatorCPU_new::get_singleton());
    h_L_nnzinsubs_splitrhs.alloc();
    if(h_L->order == 'C') {
        memcpy(h_L_nnzinsubs_splitrhs.vals, h_L->ptrs, h_L_nnzinsubs_splitrhs.size * sizeof(I));
    }
    else if(h_L_reordered.vals != nullptr && h_L_reordered.order == 'C') {
        memcpy(h_L_nnzinsubs_splitrhs.vals, h_L_reordered.ptrs, h_L_nnzinsubs_splitrhs.size * sizeof(I));
    }
    else {
        // h_L->order == 'R'
        memset(h_L_nnzinsubs_splitrhs.vals, 0, h_L_nnzinsubs_splitrhs.size * sizeof(I));
        for(size_t i = 0; i < h_L->nnz; i++) {
            I col = h_L->idxs[i];
            h_L_nnzinsubs_splitrhs.vals[col]++;
        }
        I curr = 0;
        for(size_t i = 0; i < h_L_nnzinsubs_splitrhs.size; i++) {
            I val = h_L_nnzinsubs_splitrhs.vals[i];
            h_L_nnzinsubs_splitrhs.vals[i] = curr;
            curr += val;
        }
    }
    for(size_t i = 0; i < h_L_nnzinsubs_splitrhs.size; i++) {
        h_L_nnzinsubs_splitrhs.vals[i] = h_L->nnz - h_L_nnzinsubs_splitrhs.vals[i];
    }

    d_L_sp.set(h_L_to_use->nrows, h_L_to_use->ncols, h_L_to_use->nnz, h_L_to_use->order, ator_ws_tmp_linear.get());
    d_L_sp.prop.uplo = h_L_to_use->prop.uplo;
    d_L_sp.prop.diag = h_L_to_use->prop.diag;
    wss_tmp_preprocess_linear += d_L_sp.get_memory_impact();
    wss_tmp_perform_linear += d_L_sp.get_memory_impact();

    ops_chunks.resize(num_chunks);

    first_dense_chunk = num_chunks;

    for(size_t ch = 0; ch < num_chunks; ch++) {
        stacktimer::info("trsm_hcsx_ddny_tri_splitrhs::setup chunk %zu", ch);

        size_t rhs_start = partition.vals[ch];
        size_t rhs_end = partition.vals[ch+1];

        char spdn_factor = '_';
        if(cfg.spdn_criteria == 'S') {
            spdn_factor = 'S';
        }
        if(cfg.spdn_criteria == 'D') {
            spdn_factor = 'D';
        }
        if(cfg.spdn_criteria == 'C') {
            double fraction_treshold = cfg.spdn_param;
            double curr_fraction = ((double)ch + 0.5) / num_chunks;
            spdn_factor = ((curr_fraction < fraction_treshold) ? 'S' : 'D');
        }
        if(cfg.spdn_criteria == 'Z') {
            double fraction_treshold = cfg.spdn_param;
            double curr_fraction = (double)(rhs_start + rhs_end) / 2.0 / d_X->ncols;
            spdn_factor = ((curr_fraction < fraction_treshold) ? 'S' : 'D');
        }
        if(cfg.spdn_criteria == 'T') {
            size_t k_start = h_X_colpivots.vals[rhs_start];
            size_t k_size = h_L_to_use->nrows - k_start;
            math::operations::submatrix_csx_csy<T,I> op_sub_L;
            op_sub_L.set_matrix_src(h_L_to_use);
            op_sub_L.set_bounds(k_start, h_L_to_use->nrows, k_start, h_L_to_use->ncols);
            op_sub_L.setup();
            size_t nnz = op_sub_L.get_output_matrix_nnz() + (h_L_to_use->prop.diag == 'U') * k_size;
            size_t nvals = k_size * (k_size + 1) / 2;
            double fraction_treshold = cfg.spdn_param;
            double curr_fraction = (double)nnz / nvals;
            spdn_factor = ((curr_fraction < fraction_treshold) ? 'S' : 'D');
        }

        if(spdn_factor == 'D' && first_dense_chunk == num_chunks) {
            first_dense_chunk = ch;
            
            size_t dense_rhs_start = partition.vals[first_dense_chunk];
            size_t dense_k_start = h_X_colpivots.vals[dense_rhs_start];
            size_t dense_k_size = h_L_to_use->nrows - dense_k_start;
            d_L_dn.set(dense_k_size, dense_k_size, cfg.factor_order_dn, ator_ws_tmp_linear.get());
            d_L_dn.prop.uplo = d_L_sp.prop.uplo;
            d_L_dn.prop.diag = d_L_sp.prop.diag;
            wss_tmp_perform_linear += d_L_dn.get_memory_impact();

            op_d_sub_L_sp2dn = submatrix_dcsx_ddny<T,I>::make();
            op_d_sub_L_sp2dn->set_handles(q);
            op_d_sub_L_sp2dn->set_bounds(dense_k_start, d_L_sp.nrows, dense_k_start, d_L_sp.ncols);
            op_d_sub_L_sp2dn->set_matrix_src(&d_L_sp);
            op_d_sub_L_sp2dn->set_matrix_dst(&d_L_dn);
            op_d_sub_L_sp2dn->setup();
            wss_internal += op_d_sub_L_sp2dn->get_wss_internal();
            wss_persistent += utils::round_up(op_d_sub_L_sp2dn->get_wss_persistent(), ator_ws_persistent->get_align());
            wss_tmp_preprocess_overlap = std::max(wss_tmp_preprocess_overlap, op_d_sub_L_sp2dn->get_wss_tmp_preprocess());
            wss_tmp_perform_overlap = std::max(wss_tmp_perform_overlap, op_d_sub_L_sp2dn->get_wss_tmp_perform());
        }

        gpu_trsm_trirhs_chunk_splitrhs<T,I> & op_chunk = ops_chunks[ch];
        op_chunk.set_config(spdn_factor);
        op_chunk.set_range(rhs_start, rhs_end);
        op_chunk.set_handles(q, handle_spblas, handle_dnblas);
        op_chunk.set_matrix_d_L_sp(&d_L_sp);
        op_chunk.set_matrix_d_L_dn(&d_L_dn);
        op_chunk.set_matrix_d_X(d_X);
        op_chunk.set_h_X_colpivots(&h_X_colpivots);
        op_chunk.set_h_L_nnzinsubs(&h_L_nnzinsubs_splitrhs);
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

    stacktimer::pop();

    called_setup = true;
}



template<typename T, typename I>
size_t trsm_hcsx_ddny_tri_splitrhs<T,I>::get_wss_internal()
{
    if(!called_setup) eslog::error("setup was not called\n");

    return wss_internal;
}



template<typename T, typename I>
size_t trsm_hcsx_ddny_tri_splitrhs<T,I>::get_wss_persistent()
{
    if(!called_setup) eslog::error("setup was not called\n");

    return wss_persistent;
}



template<typename T, typename I>
size_t trsm_hcsx_ddny_tri_splitrhs<T,I>::get_wss_tmp_preprocess()
{
    if(!called_setup) eslog::error("setup was not called\n");

    return wss_tmp_preprocess;
}



template<typename T, typename I>
size_t trsm_hcsx_ddny_tri_splitrhs<T,I>::get_wss_tmp_perform()
{
    if(!called_setup) eslog::error("setup was not called\n");

    return wss_tmp_perform;
}



template<typename T, typename I>
void trsm_hcsx_ddny_tri_splitrhs<T,I>::set_ws_persistent(void * ws_persistent_)
{
    if(ws_persistent_ == nullptr && wss_persistent > 0) eslog::error("persistent workspace is null\n");
    if(ws_persistent != nullptr) eslog::error("cannot re-set persistent workspace\n");

    ws_persistent = ws_persistent_;
}



template<typename T, typename I>
void trsm_hcsx_ddny_tri_splitrhs<T,I>::preprocess_submit(void * ws_tmp)
{
    stacktimer::push("trsm_hcsx_ddny_tri_splitrhs::preprocess_submit");

    if(!called_setup) eslog::error("setup has not been called\n");
    if(called_preprocess) eslog::error("preprocess has already been called\n");
    if(ws_tmp == nullptr && wss_tmp_preprocess > 0) eslog::error("temporary workspace is null\n");

    ator_ws_persistent->set(ws_persistent, wss_persistent);

    ator_ws_tmp_linear->set(ws_tmp, wss_tmp_preprocess_linear);
    ator_ws_tmp_overlap->set((char*)ws_tmp + wss_tmp_preprocess_linear, wss_tmp_preprocess_overlap);

    d_L_sp.alloc();
    gpu::mgm::copy_submit(q, *h_L_to_use, d_L_sp, true, true);

    for(size_t ch = 0; ch < num_chunks; ch++) {
        stacktimer::info("trsm_hcsx_ddny_tri_splitrhs::preprocess_submit chunk %zu", ch);

        gpu_trsm_trirhs_chunk_splitrhs<T,I> & op_chunk = ops_chunks[ch];
        op_chunk.set_ws_persistent(ator_ws_persistent->alloc(op_chunk.get_wss_persistent()));
        op_chunk.preprocess_submit(ator_ws_tmp_overlap->alloc(op_chunk.get_wss_tmp_preprocess()));
    }

    if(first_dense_chunk < num_chunks) {
        op_d_sub_L_sp2dn->set_ws_persistent(ator_ws_persistent->alloc(op_d_sub_L_sp2dn->get_wss_persistent()));
        op_d_sub_L_sp2dn->preprocess_submit(ator_ws_tmp_overlap->alloc(op_d_sub_L_sp2dn->get_wss_tmp_preprocess()));
    }

    d_L_sp.free();
    
    ator_ws_tmp_linear->unset();
    ator_ws_tmp_overlap->unset();

    stacktimer::pop();

    called_preprocess = true;
}



template<typename T, typename I>
void trsm_hcsx_ddny_tri_splitrhs<T,I>::perform_submit(void * ws_tmp)
{
    stacktimer::push("trsm_hcsx_ddny_tri_splitrhs::perform_submit");

    if(!called_preprocess) eslog::error("preprocess has not been called\n");
    if(ws_tmp == nullptr && wss_tmp_perform > 0) eslog::error("temporary workspace is null\n");

    ator_ws_tmp_linear->set(ws_tmp, wss_tmp_perform_linear);
    ator_ws_tmp_overlap->set((char*)ws_tmp + wss_tmp_perform_linear, wss_tmp_perform_overlap);

    if(cfg.factor_order_sp != h_L->order || cfg.factor_order_dn != h_L->order) {
        op_h_L_reorder.perform_values();
    }

    d_L_sp.alloc();
    gpu::mgm::copy_submit(q, *h_L_to_use, d_L_sp, true, true);

    if(first_dense_chunk < num_chunks) {
        d_L_dn.alloc();
        op_d_sub_L_sp2dn->perform_submit(ator_ws_tmp_overlap->alloc(op_d_sub_L_sp2dn->get_wss_tmp_perform()));
    }

    for(size_t ch = 0; ch < num_chunks; ch++) {
        stacktimer::info("trsm_hcsx_ddny_tri_splitrhs::perform_submit chunk %zu", ch);

        gpu_trsm_trirhs_chunk_splitrhs<T,I> & op_chunk = ops_chunks[ch];
        op_chunk.perform_submit(ator_ws_tmp_overlap->alloc(op_chunk.get_wss_tmp_perform()));
    }

    d_L_dn.free();
    d_L_sp.free();

    ator_ws_tmp_linear->unset();
    ator_ws_tmp_overlap->unset();

    stacktimer::pop();
}



#define INSTANTIATE_T_I(T,I) \
template class trsm_hcsx_ddny_tri_splitrhs<T,I>;

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
