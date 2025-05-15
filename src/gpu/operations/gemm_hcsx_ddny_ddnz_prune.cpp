
#include "gpu/operations/gemm_hcsx_ddny_ddnz_prune.h"

#include "basis/utilities/stacktimer.h"



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
void gemm_hcsx_ddny_ddnz_prune<T,I>::set_handles(gpu::mgm::queue q_, gpu::spblas::handle handle_spblas_, gpu::dnblas::handle handle_dnblas_)
{
    if(called_set_handles) eslog::error("handles are already set\n");

    q = q_;
    handle_spblas = handle_spblas_;
    handle_dnblas = handle_dnblas_;

    called_set_handles = true;
}



template<typename T, typename I>
void gemm_hcsx_ddny_ddnz_prune<T,I>::set_matrix_h_A(MatrixCsxView_new<T,I> * h_A_)
{
    if(h_A != nullptr) eslog::error("matrix h_A is already set\n");
    if(h_A_ == nullptr) eslog::error("h_A cannot be nullptr\n");

    h_A = h_A_;
}



template<typename T, typename I>
void gemm_hcsx_ddny_ddnz_prune<T,I>::set_matrix_d_B(MatrixDenseView_new<T> * d_B_)
{
    if(d_B != nullptr) eslog::error("matrix d_B is already set\n");
    if(d_B_ == nullptr) eslog::error("d_B cannot be nullptr\n");

    d_B = d_B_;
}



template<typename T, typename I>
void gemm_hcsx_ddny_ddnz_prune<T,I>::set_matrix_d_C(MatrixDenseView_new<T> * d_C_)
{
    if(d_C != nullptr) eslog::error("matrix d_C is already set\n");
    if(d_C_ == nullptr) eslog::error("d_C cannot be nullptr\n");

    d_C = d_C_;
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
    stacktimer::push("gemm_hcsx_ddny_ddnz_prune::setup");

    if(!called_set_config) eslog::error("config is not set\n");
    if(!called_set_handles) eslog::error("handles are not set\n");
    if(h_A == nullptr) eslog::error("matrix A is not set\n");
    if(d_B == nullptr) eslog::error("matrix B is not set\n");
    if(d_C == nullptr) eslog::error("matrix C is not set\n");
    if(!h_A->ator->is_data_accessible_cpu()) eslog::error("matrix h_A must be cpu-accessible\n");
    if(!d_B->ator->is_data_accessible_gpu()) eslog::error("matrix d_B must be gpu-accessible\n");
    if(!d_C->ator->is_data_accessible_gpu()) eslog::error("matrix d_C must be gpu-accessible\n");
    if(called_setup) eslog::error("setup has already been called\n");
    if(h_A->nrows != d_C->nrows || d_B->ncols != d_C->ncols || h_A->ncols != d_B->nrows) eslog::error("incompatible matrices\n");

    ator_persistent = std::make_unique<AllocatorArena_new>(AllocatorGPU_new::get_singleton());
    ator_tmp_linear = std::make_unique<AllocatorArena_new>(AllocatorGPU_new::get_singleton());
    ator_tmp_overlap = std::make_unique<AllocatorSinglePointer_new>(AllocatorGPU_new::get_singleton());

    if(prune_rows || prune_cols) {
        op_h_prune_A.set_matrix_src(h_A);
        op_h_prune_A.set_matrix_dst_sp(&h_A_pruned);
        op_h_prune_A.set_pruning_mode(prune_rows, prune_cols);
        op_h_prune_A.setup();

        m = op_h_prune_A.get_dst_matrix_nrows();
        n = d_C->ncols;
        k = op_h_prune_A.get_dst_matrix_ncols();

        h_A_pruned.set(m, k, h_A->nnz, h_A->order, AllocatorHostPinned_new::get_singleton());

        if(prune_rows) {
            h_pruned_rows.set(m, AllocatorHostPinned_new::get_singleton());
            op_h_prune_A.set_vector_pruned_rows(&h_pruned_rows);
        }
        if(prune_cols) {
            h_pruned_cols.set(k, AllocatorHostPinned_new::get_singleton());
            op_h_prune_A.set_vector_pruned_cols(&h_pruned_cols);
        }

        h_A_to_use = &h_A_pruned;
    }
    else {
        m = d_C->nrows;
        n = d_C->ncols;
        k = h_A->ncols;

        h_A_to_use = h_A;
    }

    if(prune_rows) {
        d_pruned_rows.set(h_pruned_rows.size, ator_persistent.get());
        wss_persistent += utils::round_up(d_pruned_rows.get_memory_impact(), ator_persistent->get_align());
    }
    if(prune_cols) {
        d_pruned_cols.set(h_pruned_cols.size, ator_persistent.get());
        wss_persistent += utils::round_up(d_pruned_cols.get_memory_impact(), ator_persistent->get_align());
    }

    d_A_pruned_sp.set(m, k, h_A_to_use->nnz, h_A_to_use->order, ator_tmp_linear.get());
    wss_tmp_preprocess_linear += d_A_pruned_sp.get_memory_impact();
    wss_tmp_perform_linear += d_A_pruned_sp.get_memory_impact();

    if(spdn_A == 'D') {
        d_A_pruned_dn.set(m, k, h_A_to_use->order, ator_tmp_linear.get());
        wss_tmp_perform_linear += d_A_pruned_dn.get_memory_impact();
    }

    if(prune_cols) {
        d_B_pruned.set(k, n, d_B->order, ator_tmp_linear.get());
        wss_tmp_perform_linear += d_B_pruned.get_memory_impact();
        d_B_to_use = &d_B_pruned;
    }
    else {
        d_B_to_use = d_B;
    }

    if(prune_rows) {
        d_C_pruned.set(m, n, d_C->order, ator_tmp_linear.get());
        wss_tmp_perform_linear += d_C_pruned.get_memory_impact();
        d_C_to_use = &d_C_pruned;
    }
    else {
        d_C_to_use = d_C;
    }

    if(spdn_A == 'D') {
        op_d_sp2dn_A = convert_dcsx_ddny<T,I>::make();
        op_d_sp2dn_A->set_handles(q, handle_spblas);
        op_d_sp2dn_A->set_matrix_src(&d_A_pruned_sp);
        op_d_sp2dn_A->set_matrix_dst(&d_A_pruned_dn);
        op_d_sp2dn_A->setup();
        wss_internal += op_d_sp2dn_A->get_wss_internal();
        wss_persistent += utils::round_up(op_d_sp2dn_A->get_wss_persistent(), ator_persistent->get_align());
        wss_tmp_preprocess_overlap = std::max(wss_tmp_preprocess_overlap, op_d_sp2dn_A->get_wss_tmp_preprocess());
        wss_tmp_perform_overlap = std::max(wss_tmp_perform_overlap, op_d_sp2dn_A->get_wss_tmp_perform());
    }

    if(prune_cols) {
        op_d_sub_B = submatrix_ddnx_ddnx_noncontig<T,I>::make();
        op_d_sub_B->set_handles(q);
        op_d_sub_B->set_matrix_src(d_B);
        op_d_sub_B->set_matrix_dst(&d_B_pruned);
        op_d_sub_B->set_row_map(&d_pruned_cols);
        op_d_sub_B->set_col_map(nullptr);
    }
    if(prune_rows) {
        op_d_sub_C = submatrix_ddnx_ddnx_noncontig<T,I>::make();
        op_d_sub_C->set_handles(q);
        op_d_sub_C->set_matrix_src(d_C);
        op_d_sub_C->set_matrix_dst(&d_C_pruned);
        op_d_sub_C->set_row_map(&d_pruned_rows);
        op_d_sub_C->set_col_map(nullptr);

        op_d_super_C = supermatrix_ddnx_ddnx_noncontig<T,I>::make();
        op_d_super_C->set_handles(q);
        op_d_super_C->set_matrix_src(&d_C_pruned);
        op_d_super_C->set_matrix_dst(d_C);
        op_d_super_C->set_row_map(&d_pruned_rows);
        op_d_super_C->set_col_map(nullptr);
    }

    if(spdn_A == 'S') {
        op_gemm_sp = gemm_dcsx_ddny_ddnz<T,I>::make();
        op_gemm_sp->set_handles(q, handle_spblas);
        op_gemm_sp->set_matrix_A(&d_A_pruned_sp);
        op_gemm_sp->set_matrix_B(d_B_to_use);
        op_gemm_sp->set_matrix_C(d_C_to_use);
        op_gemm_sp->set_coefficients(alpha, beta);
        op_gemm_sp->setup();
        wss_internal += op_gemm_sp->get_wss_internal();
        wss_persistent += utils::round_up(op_gemm_sp->get_wss_persistent(), ator_persistent->get_align());
        wss_tmp_preprocess_overlap = std::max(wss_tmp_preprocess_overlap, op_gemm_sp->get_wss_tmp_preprocess());
        wss_tmp_perform_overlap = std::max(wss_tmp_perform_overlap, op_gemm_sp->get_wss_tmp_perform());
    }
    if(spdn_A == 'D') {
        op_gemm_dn = gemm_ddnx_ddny_ddnz<T>::make();
        op_gemm_dn->set_handles(q, handle_dnblas);
        op_gemm_dn->set_matrix_A(&d_A_pruned_dn);
        op_gemm_dn->set_matrix_B(d_B_to_use);
        op_gemm_dn->set_matrix_C(d_C_to_use);
        op_gemm_dn->set_coefficients(alpha, beta);
        op_gemm_dn->setup();
        wss_tmp_perform_overlap = std::max(wss_tmp_perform_overlap, op_gemm_dn->get_wss_tmp_perform());
    }

    wss_tmp_preprocess_linear = utils::round_up(wss_tmp_preprocess_linear, ator_tmp_linear->get_align());
    wss_tmp_perform_linear = utils::round_up(wss_tmp_perform_linear, ator_tmp_overlap->get_align());

    wss_tmp_preprocess = wss_tmp_preprocess_linear + wss_tmp_preprocess_overlap;
    wss_tmp_perform = wss_tmp_perform_linear + wss_tmp_perform_overlap;

    // stacktimer::info("wss_internal       %zu", wss_internal);
    // stacktimer::info("wss_persistent     %zu", wss_persistent);
    // stacktimer::info("wss_tmp_preprocess %zu", wss_tmp_preprocess);
    // stacktimer::info("wss_tmp_perform    %zu", wss_tmp_perform);

    stacktimer::pop();

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
    stacktimer::push("gemm_hcsx_ddny_ddnz_prune::preprocess_submit");

    if(!called_setup) eslog::error("setup has not been called\n");
    if(called_preprocess) eslog::error("preprocess has already been called\n");
    if(ws_tmp == nullptr && wss_tmp_preprocess > 0) eslog::error("temporary workspace is null\n");

    ator_persistent->set(ws_persistent, wss_persistent);

    ator_tmp_linear->set(ws_tmp, wss_tmp_preprocess_linear);
    ator_tmp_overlap->set((char*)ws_tmp + wss_tmp_preprocess_linear, wss_tmp_preprocess_overlap);

    if(prune_rows) {
        h_pruned_rows.alloc();
    }
    if(prune_cols) {
        h_pruned_cols.alloc();
    }

    if(prune_rows || prune_cols) {
        h_A_pruned.alloc();
        op_h_prune_A.preprocess();
        op_h_prune_A.perform();
    }

    if(prune_rows) {
        d_pruned_rows.alloc();
        gpu::mgm::copy_submit(q, h_pruned_rows, d_pruned_rows);
    }
    if(prune_cols) {
        d_pruned_cols.alloc();
        gpu::mgm::copy_submit(q, h_pruned_cols, d_pruned_cols);
    }

    d_A_pruned_sp.alloc();
    gpu::mgm::copy_submit(q, *h_A_to_use, d_A_pruned_sp, true, false);

    if(spdn_A == 'D') {
        op_d_sp2dn_A->set_ws_persistent(ator_persistent->alloc(op_d_sp2dn_A->get_wss_persistent()));
        op_d_sp2dn_A->preprocess_submit(ator_tmp_overlap->alloc(op_d_sp2dn_A->get_wss_tmp_preprocess()));
    }

    if(spdn_A == 'S') {
        op_gemm_sp->set_ws_persistent(ator_persistent->alloc(op_gemm_sp->get_wss_persistent()));
        op_gemm_sp->preprocess_submit(ator_tmp_overlap->alloc(op_gemm_sp->get_wss_tmp_preprocess()));
    }

    d_A_pruned_sp.free();

    ator_tmp_linear->unset();
    ator_tmp_overlap->unset();

    stacktimer::pop();

    called_preprocess = true;
}



template<typename T, typename I>
void gemm_hcsx_ddny_ddnz_prune<T,I>::perform_submit(void * ws_tmp)
{
    stacktimer::push("gemm_hcsx_ddny_ddnz_prune::perform_submit");

    if(!called_preprocess) eslog::error("preprocess has not been called\n");
    if(ws_tmp == nullptr && wss_tmp_perform > 0) eslog::error("temporary workspace is null\n");

    ator_tmp_linear->set(ws_tmp, wss_tmp_perform_linear);
    ator_tmp_overlap->set((char*)ws_tmp + wss_tmp_perform_linear, wss_tmp_perform_overlap);

    d_A_pruned_sp.alloc();
    if(spdn_A == 'D') d_A_pruned_dn.alloc();
    if(d_B_to_use == &d_B_pruned) d_B_pruned.alloc();
    if(d_C_to_use == &d_C_pruned) d_C_pruned.alloc();

    if(prune_rows || prune_cols) {
        op_h_prune_A.perform();
    }

    gpu::mgm::copy_submit(q, *h_A_to_use, d_A_pruned_sp, true, true);

    if(spdn_A == 'D') {
        op_d_sp2dn_A->perform_submit(ator_tmp_overlap->alloc(op_d_sp2dn_A->get_wss_tmp_perform()));
    }

    if(prune_cols) {
        op_d_sub_B->perform_submit();
    }
    if(prune_rows) {
        op_d_sub_C->perform_submit();
    }

    if(spdn_A == 'S') {
        op_gemm_sp->perform_submit(ator_tmp_overlap->alloc(op_gemm_sp->get_wss_tmp_perform()));
    }
    if(spdn_A == 'D') {
        op_gemm_dn->perform_submit(ator_tmp_overlap->alloc(op_gemm_dn->get_wss_tmp_perform()));
    }

    if(prune_rows) {
        op_d_super_C->perform_submit();
    }

    d_A_pruned_sp.free();
    if(spdn_A == 'D') d_A_pruned_dn.free();
    if(d_B_to_use == &d_B_pruned) d_B_pruned.free();
    if(d_C_to_use == &d_C_pruned) d_C_pruned.free();

    ator_tmp_linear->unset();
    ator_tmp_overlap->unset();

    stacktimer::pop();
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
        INSTANTIATE_T(std::complex<double>)

            INSTANTIATE

        #undef INSTANTIATE
    #undef INSTANTIATE_T
#undef INSTANTIATE_T_I



}
}
}
