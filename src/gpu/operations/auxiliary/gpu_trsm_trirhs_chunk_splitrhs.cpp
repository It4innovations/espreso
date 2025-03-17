
#include "gpu/operations/auxiliary/gpu_trsm_trirhs_chunk_splitrhs.h"



namespace espreso {
namespace gpu {
namespace operations {



template<typename T, typename I>
void gpu_trsm_trirhs_chunk_splitrhs<T,I>::set_config(config cfg_)
{
    if(called_set_config) eslog::error("config has already been set\n");

    factor_spdn = factor_spdn_;

    called_set_config = true;
}



template<typename T, typename I>
void gpu_trsm_trirhs_chunk_splitrhs<T,I>::set_range(size_t rhs_start_, size_t rhs_end_)
{
    if(called_set_range) eslog::error("range has already been set\n");

    rhs_start = rhs_start_;
    rhs_end = rhs_end_;
    rhs_size = rhs_size_;

    called_set_range = true;
}



template<typename T, typename I>
void gpu_trsm_trirhs_chunk_splitrhs<T,I>::set_handles(gpu::mgm::queue q_, gpu::spblas::handle handle_spblas_, gpu::dnblas::handle handle_dnblas_)
{
    if(called_set_handles) eslog::error("handles have already been set\n");

    q = q_;
    handle_spblas = handle_spblas_;
    handle_dnblas = handle_dnblas_;

    called_set_handles = true;
}



template<typename T, typename I>
void gpu_trsm_trirhs_chunk_splitrhs<T,I>::set_matrix_d_L_sp(MatrixCsxView_new<T,I> * d_L_sp_)
{
    if(d_L_sp != nullptr) eslog::error("matrix d_L_sp is already set\n");
    if(d_L_sp_ == nullptr) eslog::error("d_L_sp cannot be nullptr\n");

    d_L_sp = d_L_sp_;
}



template<typename T, typename I>
void gpu_trsm_trirhs_chunk_splitrhs<T,I>::set_matrix_d_L_dn(MatrixDenseView_new<T> * d_L_dn_)
{
    if(d_L_dn != nullptr) eslog::error("matrix d_L_dn is already set\n");
    if(d_L_dn_ == nullptr) eslog::error("d_L_dn cannot be nullptr\n");

    d_L_dn = d_L_dn_;
}



template<typename T, typename I>
void gpu_trsm_trirhs_chunk_splitrhs<T,I>::set_matrix_d_X(MatrixDenseView_new<T> * d_X_)
{
    if(d_X != nullptr) eslog::error("matrix d_X is already set\n");
    if(d_X_ == nullptr) eslog::error("d_X cannot be nullptr\n");

    d_X = d_X_;
}



template<typename T, typename I>
void gpu_trsm_trirhs_chunk_splitrhs<T,I>::set_h_X_colpivots(VectorDenseView_new<I> * h_X_colpivots_)
{
    if(h_X_colpivots != nullptr) eslog::error("X colpivots are already set\n");
    if(h_X_colpivots_ == nullptr) eslog::error("X colpivots cannot be nullptr\n");

    h_X_colpivots = h_X_colpivots_;
}



template<typename T, typename I>
void gpu_trsm_trirhs_chunk_splitrhs<T,I>::set_h_L_nnzinsubs(VectorDenseView_new<I> * h_L_nnzinsubs_)
{
    if(h_L_nnzinsubs != nullptr) eslog::error("Lnnzinsubs is already set\n");
    if(h_L_nnzinsubs_ == nullptr) eslog::error("Lnnzinsubs cannot be nullptr\n");

    h_L_nnzinsubs = h_L_nnzinsubs_;
}



template<typename T, typename I>
void gpu_trsm_trirhs_chunk_splitrhs<T,I>::setup()
{
    if(!called_set_config) eslog::error("config is not set\n");
    if(!called_set_range) eslog::error("range is not set\n");
    if(!called_set_handles) eslog::error("handles are not set\n");
    if(d_L_sp == nullptr) eslog::error("matrix d_L_sp is not set\n");
    if(d_L_dn == nullptr) eslog::error("matrix d_L_dn is not set\n");
    if(d_X == nullptr) eslog::error("matrix X is not set\n");
    if(h_X_colpivots == nullptr) eslog::error("X colpivots are not set\n");
    if(h_L_nnzinsubs == nullptr) eslog::error("Lnnzinsubs are not set\n");
    if(d_L_sp->nrows != d_L_sp->ncols) eslog::error("Lsp is not square\n");
    if(d_L_dn->nrows != d_L_dn->ncols) eslog::error("Ldn is not square\n");
    if(d_L_sp->nrows != X->nrows) eslog::error("incompatible matrices\n");
    if(h_X_colpivots->size != X->nrows) eslog::error("wrong X rowtrails size\n");
    if(h_L_nnzinsubs->size != d_L_sp->ncols) eslog::error("wrong size of Lnnzinsubs\n");
    if(d_L_sp->prop.uplo != 'L') eslog::error("matrix L must have uplo=L\n");

    ator_ws_persistent = std::make_unique<AllocatorArena_new>(false, true, gpu::mgm::get_natural_pitch_align());
    ator_ws_tmp_linear = std::make_unique<AllocatorArena_new>(false, true, gpu::mgm::get_natural_pitch_align());
    ator_ws_tmp_overlap = std::make_unique<AllocatorSinglePointer_new>(false, true, gpu::mgm::get_natural_pitch_align());

    k_start = h_X_colpivots[rhs_start];
    k_size = d_X->nrows - k_start;

    if(cfg.factor_spdn == 'S') {
        d_sub_L_sp.set(k_size, k_size, h_L_nnzinsubs[k_start], d_L_sp->order, ator_ws_tmp_linear.get());
        d_sub_L_sp.prop.uplo = d_L_sp->prop.uplo;
        d_sub_L_sp.prop.diag = d_L_sp->prop.diag;
        wss_tmp_perform_linear += d_sub_L_sp.get_memory_impact();
    }
    if(cfg.factor_spdn == 'D') {
        d_sub_L_dn.set_view(k_size, k_size, d_L_dn->order, nullptr);
    }

    if(cfg.factor_spdn == 'S') {
        op_d_sub_L_sp.set_handles(q);
        op_d_sub_L_sp.set_bounds(k_start, d_L->nrows, k_start, d_L->ncols);
        op_d_sub_L_sp.set_matrix_src(d_L);
        op_d_sub_L_sp.set_matrix_dst(d_sub_L_sp);
        op_d_sub_L_sp.setup();
        wss_internal += op_d_sub_L_sp.get_wss_internal();
        wss_persistent += op_d_sub_L_sp.get_wss_persistent();
        wss_tmp_preprocess_overlap = std::max(wss_tmp_preprocess_overlap, op_d_sub_L_sp.get_wss_tmp_preprocess());
        wss_tmp_peform_overlap = std::max(wss_tmp_peform_overlap, op_d_sub_L_sp.get_wss_tmp_perform());
    }
    if(cfg.factor_spdn == 'D') {
        size_t dense_size = d_L_dn->nrows;
        size_t dense_start = dense_size - k_size;
        op_sub_L_dn.set_matrix_src(d_L_dn);
        op_sub_L_dn.set_matrix_dst(d_sub_L_dn);
        op_sub_L_dn.set_bounds(dense_start, d_L_dn->nrows, dense_start, d_L_dn->ncols);
    }

    d_sub_X.set_view(k_size, rhs_size, d_X->ld, d_X->order, nullptr);

    op_sub_X.set_matrix_src(d_X);
    op_sub_X.set_matrix_dst(d_sub_X);
    op_sub_X.set_bounds(k_start, d_X->nrows, rhs_start, rhs_end);

    if(cfg.factor_spdn == 'S') {
        op_d_trsm_sp.set_handles(q, handle_spblas);
        op_d_trsm_sp.set_matrix_A(d_sub_L_sp);
        op_d_trsm_sp.set_matrix_X(d_sub_X);
        op_d_trsm_sp.set_matrix_B(d_sub_X);
        op_d_trsm_sp.setup();
        wss_internal += op_d_trsm_sp.get_wss_internal();
        wss_persistent += op_d_trsm_sp.get_wss_persistent();
        wss_tmp_preprocess_overlap = std::max(wss_tmp_preprocess_overlap, op_d_trsm_sp.get_wss_tmp_preprocess());
        wss_tmp_peform_overlap = std::max(wss_tmp_peform_overlap, op_d_trsm_sp.get_wss_tmp_perform());
    }
    if(cfg.factor_spdn == 'D') {
        op_d_trsm_dn.set_handles(q, handle_dnblas);
        op_d_trsm_dn.set_matrix_A(d_sub_L_dn);
        op_d_trsm_dn.set_matrix_X(d_sub_X);
        op_d_trsm_dn.setup();
        wss_tmp_peform_overlap = std::max(wss_tmp_peform_overlap, op_d_trsm_dn.get_wss_tmp_perform());
    }

    wss_tmp_preprocess_linear = ((wss_tmp_preprocess_linear - 1) / ator_ws_tmp_linear->get_align() + 1) * ator_ws_tmp_linear->get_align();
    wss_tmp_perform_linear = ((wss_tmp_perform_linear - 1) / ator_ws_tmp_linear->get_align() + 1) * ator_ws_tmp_linear->get_align();

    wss_tmp_preprocess = wss_tmp_preprocess_linear + wss_tmp_preprocess_overlap;
    wss_tmp_perform = wss_tmp_perform_linear + wss_tmp_perform_overlap;

    called_setup = true;
}



template<typename T, typename I>
size_t gpu_trsm_trirhs_chunk_splitrhs<T,I>::get_wss_internal()
{
    if(!called_setup) eslog::error("setup was not called\n");

    return wss_internal;
}



template<typename T, typename I>
size_t gpu_trsm_trirhs_chunk_splitrhs<T,I>::get_wss_persistent()
{
    if(!called_setup) eslog::error("setup was not called\n");

    return wss_persistent;
}



template<typename T, typename I>
size_t gpu_trsm_trirhs_chunk_splitrhs<T,I>::get_wss_tmp_preprocess()
{
    if(!called_setup) eslog::error("setup was not called\n");

    return wss_tmp_preprocess;
}



template<typename T, typename I>
size_t gpu_trsm_trirhs_chunk_splitrhs<T,I>::get_wss_tmp_perform()
{
    if(!called_setup) eslog::error("setup was not called\n");

    return wss_tmp_perform;
}



template<typename T, typename I>
void gpu_trsm_trirhs_chunk_splitrhs<T,I>::set_ws_persistent(void * ws_persistent_)
{
    if(ws_persistent_ == nullptr && wss_persistent > 0) eslog::error("persistent workspace is null\n");
    if(ws_persistent != nullptr) eslog::error("cannot re-set persistent workspace\n");

    ws_persistent = ws_persistent_;
}



template<typename T, typename I>
void gpu_trsm_trirhs_chunk_splitrhs<T,I>::preprocess_submit(void * ws_tmp)
{
    if(!called_setup) eslog::error("setup has not been called\n");
    if(called_preprocess) eslog::error("preprocess has already been called\n");
    if(ws_tmp == nullptr && wss_tmp_preprocess > 0) eslog::error("temporary workspace is null\n");

    ator_ws_persistent.set(ws_persistent, wss_persistent);

    ator_ws_tmp_linear.set(ws_tmp, wss_tmp_preprocess_linear);
    ator_ws_tmp_overlap.set(ws_tmp + wss_tmp_preprocess_linear, wss_tmp_preprocess_overlap);

    if(cfg.factor_spdn == 'S') {
        op_d_sub_L_sp.set_ws_persistent(ator_ws_persistent.alloc(op_d_sub_L_sp.get_wss_persistent()));
        op_d_sub_L_sp.preprocess_submit(ator_ws_tmp_overlap.alloc(op_d_sub_L_sp.get_wss_tmp_preprocess()));
    }

    if(cfg.factor_spdn == 'S') {
        op_d_trsm_sp.set_ws_persistent(ator_ws_persistent.alloc(op_d_trsm_sp.get_wss_persistent()));
        op_d_trsm_sp.preprocess_submit(ator_ws_tmp_overlap.alloc(op_d_trsm_sp.get_wss_tmp_preprocess()));
    }

    ator_ws_tmp_linear.unset();
    ator_ws_tmp_overlap.unset();

    called_preprocess = true;
}



template<typename T, typename I>
void gpu_trsm_trirhs_chunk_splitrhs<T,I>::perform_submit(void * ws_tmp)
{
    if(!called_preprocess) eslog::error("preprocess has not been called\n");
    if(ws_tmp == nullptr && wss_tmp_perform > 0) eslog::error("temporary workspace is null\n");

    ator_ws_tmp_linear.set(ws_tmp, wss_tmp_perform_linear);
    ator_ws_tmp_overlap.set(ws_tmp + wss_tmp_perform_linear, wss_tmp_perform_overlap);

    if(cfg.factor_spdn == 'S') {
        d_sub_L_sp.alloc();
    }

    if(cfg.factor_spdn == 'S') {
        op_d_sub_L_sp.perform_submit(ator_ws_tmp_overlap.alloc(op_d_sub_L_sp.get_wss_tmp_perform()));
    }
    if(cfg.factor_spdn == 'S') {
        op_sub_L_dn.perform();
    }

    d_sub_X.perform();

    if(cfg.factor_spdn == 'S') {
        op_d_trsm_sp.perform_submit(ator_ws_tmp_overlap.alloc(op_d_trsm_sp.get_wss_tmp_perform()));
    }
    if(cfg.factor_spdn == 'D') {
        op_d_trsm_dn.perform_submit(ator_ws_tmp_overlap.alloc(op_d_trsm_dn.get_wss_tmp_perform()));
    }

    if(cfg.factor_spdn == 'S') {
        d_sub_L_sp.free();
    }
    
    ator_ws_tmp_linear.unset();
    ator_ws_tmp_overlap.unset();
}



#define INSTANTIATE_T_I(T,I) \
template class gpu_trsm_trirhs_chunk_splitrhs<T,I>;

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
