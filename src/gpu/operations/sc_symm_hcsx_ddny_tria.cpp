
#include "gpu/operations/sc_symm_hcsx_ddny_tria.h"

#include "basis/utilities/stacktimer.h"
#include "math/operations/pivots_trails_csx.h"
#include "math/operations/sorting_permutation.h"
#include "math/operations/permute_csx_csx.h"



namespace espreso {
namespace gpu {
namespace operations {



template<typename T, typename I>
void sc_symm_hcsx_ddny_tria<T,I>::set_config(config cfg_)
{
    if(called_set_config) eslog::error("config has already been set\n");

    cfg = cfg_;

    called_set_config = true;
}



template<typename T, typename I>
void sc_symm_hcsx_ddny_tria<T,I>::set_handles(gpu::mgm::queue q_, gpu::spblas::handle handle_spblas_, gpu::dnblas::handle handle_dnblas_)
{
    if(called_set_handles) eslog::error("handles have already been set\n");

    q = q_;
    handle_spblas = handle_spblas_;
    handle_dnblas = handle_dnblas_;

    called_set_handles = true;
}



template<typename T, typename I>
void sc_symm_hcsx_ddny_tria<T,I>::set_coefficients(Treal alpha_)
{
    alpha = alpha_;
}



template<typename T, typename I>
void sc_symm_hcsx_ddny_tria<T,I>::set_h_A11_solver(DirectSparseSolver<T,I> * h_A11_solver_)
{
    if(h_A11_solver != nullptr) eslog::error("h_A11_solver is already set\n");
    if(h_A11_solver_ == nullptr) eslog::error("h_A11_solver cannot be nullptr\n");

    h_A11_solver = h_A11_solver_;
}



template<typename T, typename I>
void sc_symm_hcsx_ddny_tria<T,I>::set_h_A12(MatrixCsxView_new<T,I> * h_A12_)
{
    if(h_A12 != nullptr) eslog::error("matrix h_A12 is already set\n");
    if(h_A12_ == nullptr) eslog::error("h_A12 cannot be nullptr\n");

    h_A12 = h_A12_;
}



template<typename T, typename I>
void sc_symm_hcsx_ddny_tria<T,I>::set_d_sc(MatrixDenseView_new<T> * d_sc_)
{
    if(d_sc != nullptr) eslog::error("matrix d_sc is already set\n");
    if(d_sc_ == nullptr) eslog::error("d_sc cannot be nullptr\n");

    d_sc = d_sc_;
}



template<typename T, typename I>
void sc_symm_hcsx_ddny_tria<T,I>::setup()
{
    stacktimer::push("sc_symm_hcsx_ddny_tria::setup");

    if(!called_set_config) eslog::error("config is not set\n");
    if(!called_set_handles) eslog::error("handles are not set\n");
    if(called_setup) eslog::error("setup has already been called\n");
    if(h_A11_solver == nullptr) eslog::error("A11 solver is not set\n");
    if(h_A12 == nullptr) eslog::error("matrix A12 is not set\n");
    if(d_sc == nullptr) eslog::error("sc matrix is not set\n");
    if(!h_A12->ator->is_data_accessible_cpu()) eslog::error("matrix h_A12 must be cpu-accessible\n");
    if(!d_sc->ator->is_data_accessible_gpu()) eslog::error("matrix d_sc must be gpu-accessible\n");
    if((size_t)h_A11_solver->getMatrixSize() != h_A12->nrows) eslog::error("incompatible matrices\n");
    if(d_sc->ncols != h_A12->ncols) eslog::error("incompatible matrices\n");
    if(d_sc->prop.uplo != 'L' && d_sc->prop.uplo != 'U') eslog::error("wrong sc uplo\n");

    if(!DirectSparseSolver<T,I>::provideFactors()) eslog::error("wrong sparse solver, must provide factors\n");
    solver_factor_uplo = '_';
    if(DirectSparseSolver<T,I>::factorsSymmetry() == Solver_Factors::HERMITIAN_LOWER) solver_factor_uplo = 'L';
    if(DirectSparseSolver<T,I>::factorsSymmetry() == Solver_Factors::HERMITIAN_UPPER) solver_factor_uplo = 'U';
    if(solver_factor_uplo == '_') eslog::error("wrong sparse solver, must be symmetric\n");

    ator_ws_persistent = std::make_unique<AllocatorArena_new>(AllocatorGPU_new::get_singleton());
    ator_ws_tmp_linear = std::make_unique<AllocatorArena_new>(AllocatorGPU_new::get_singleton());
    ator_ws_tmp_overlap = std::make_unique<AllocatorSinglePointer_new>(AllocatorGPU_new::get_singleton());

    I factor_nnz = h_A11_solver->getFactorNnz();
    A11size = h_A11_solver->getMatrixSize();

    h_factor_U_row.set(A11size, A11size, factor_nnz, 'R', AllocatorHostPinned_new::get_singleton());
    h_factor_U_row.prop.uplo = 'U';
    h_factor_U_row.prop.diag = 'N';

    h_factor_L_row.set(A11size, A11size, factor_nnz, 'R', AllocatorHostPinned_new::get_singleton());
    h_factor_L_row.prop.uplo = 'L';
    h_factor_L_row.prop.diag = 'N';

    need_reorder_factor_L2U = ((solver_factor_uplo == 'L') && (cfg.order_L == 'C'));
    need_reorder_factor_U2L = ((solver_factor_uplo == 'U') && (cfg.order_L == 'R'));

    stacktimer::push("alloc_host_factors");
    if(solver_factor_uplo == 'U' || need_reorder_factor_L2U) {
        h_factor_U_row.alloc();
    }
    if(solver_factor_uplo == 'L' || need_reorder_factor_U2L) {
        h_factor_L_row.alloc();
    }
    stacktimer::pop();

    stacktimer::push("extract_factors");
    if(solver_factor_uplo == 'U') {
        Matrix_CSR<T,I,gpu::mgm::Ah> h_factor_old_U = MatrixCsxView_new<T,I>::template to_old<gpu::mgm::Ah>(h_factor_U_row);
        h_A11_solver->getFactorU(h_factor_old_U, true, false);
    }
    if(solver_factor_uplo == 'L') {
        Matrix_CSR<T,I,gpu::mgm::Ah> h_factor_old_L = MatrixCsxView_new<T,I>::template to_old<gpu::mgm::Ah>(h_factor_L_row);
        h_A11_solver->getFactorL(h_factor_old_L, true, false);
    }
    stacktimer::pop();

    h_L_row = h_factor_L_row;
    h_L_col = h_factor_U_row.get_transposed_reordered_view();
    h_U_row = h_factor_U_row;
    h_U_col = h_factor_L_row.get_transposed_reordered_view();

    if(need_reorder_factor_L2U) {
        op_L2U.set_matrix_src(&h_L_row);
        op_L2U.set_matrix_dst(&h_L_col);
        op_L2U.perform_pattern();
    }
    if(need_reorder_factor_U2L) {
        op_U2L.set_matrix_src(&h_U_row);
        op_U2L.set_matrix_dst(&h_U_col);
        op_U2L.perform_pattern();
    }

    if(cfg.order_L == 'R') h_L_to_use = &h_L_row;
    if(cfg.order_L == 'C') h_L_to_use = &h_L_col;

    h_X_sp.set(h_A12->nrows, h_A12->ncols, h_A12->nnz, h_A12->order, AllocatorHostPinned_new::get_singleton());
    h_X_sp.alloc();

    {
        MatrixCsxData_new<T,I> h_X_sp_tmp;
        h_X_sp_tmp.set(h_X_sp.nrows, h_X_sp.ncols, h_X_sp.nnz, h_X_sp.order, AllocatorCPU_new::get_singleton());
        h_X_sp_tmp.alloc();

        {
            h_perm_fillreduce.set(h_A12->nrows, AllocatorCPU_new::get_singleton());
            h_perm_fillreduce.alloc();
            Permutation<I> h_perm_fillreduce_old = PermutationView_new<I>::template to_old<cpu_allocator>(h_perm_fillreduce);
            h_A11_solver->getPermutation(h_perm_fillreduce_old);

            math::operations::permute_csx_csx<T,I>::do_all(h_A12, &h_X_sp_tmp, &h_perm_fillreduce, nullptr);
        }

        {
            VectorDenseData_new<I> colpivots;
            colpivots.set(h_A12->ncols, AllocatorCPU_new::get_singleton());
            colpivots.alloc();

            math::operations::pivots_trails_csx<T,I>::do_all(&h_X_sp_tmp, &colpivots, 'C', 'P', '_');

            h_perm_to_sort_A12_cols.set(h_A12->ncols, AllocatorCPU_new::get_singleton());
            h_perm_to_sort_A12_cols.alloc();

            h_perm_to_sort_back_sc = h_perm_to_sort_A12_cols.get_inverse_view();

            d_perm_to_sort_back_sc.set(h_perm_to_sort_back_sc.size, ator_ws_persistent.get());
            wss_persistent += utils::round_up(d_perm_to_sort_back_sc.get_memory_impact(), ator_ws_persistent->get_align());

            math::operations::sorting_permutation<I,I>::do_all(&colpivots, &h_perm_to_sort_A12_cols);
        }

        math::operations::permute_csx_csx<T,I>::do_all(&h_X_sp_tmp, &h_X_sp, nullptr, &h_perm_to_sort_A12_cols);
    }

    d_X_sp.set(h_X_sp.nrows, h_X_sp.ncols, h_X_sp.nnz, h_X_sp.order, ator_ws_tmp_linear.get());
    wss_tmp_perform_linear += d_X_sp.get_memory_impact();

    d_X_dn.set(h_X_sp.nrows, h_X_sp.ncols, cfg.order_X, ator_ws_tmp_linear.get());
    wss_tmp_perform_linear += d_X_dn.get_memory_impact();

    d_sc_tmp1_x.set(d_sc->nrows, d_sc->ncols, d_sc->order, ator_ws_tmp_linear.get());
    d_sc_tmp1_x.prop.uplo = d_sc->prop.uplo;
    wss_tmp_perform_linear += d_sc_tmp1_x.get_memory_impact();
    d_sc_tmp1_y = d_sc_tmp1_x.get_transposed_reordered_view();

    d_sc_tmp2_x.set(d_sc->nrows, d_sc->ncols, d_sc->order, ator_ws_tmp_linear.get());
    d_sc_tmp2_x.prop.uplo = d_sc->prop.uplo;
    wss_tmp_perform_linear += d_sc_tmp2_x.get_memory_impact();
    d_sc_tmp2_y = d_sc_tmp2_x.get_transposed_reordered_view();

    op_X_sp2dn = convert_dcsx_ddny<T,I>::make();
    op_X_sp2dn->set_handles(q, handle_spblas);
    op_X_sp2dn->set_matrix_src(&d_X_sp);
    op_X_sp2dn->set_matrix_dst(&d_X_dn);
    op_X_sp2dn->setup();
    wss_internal += op_X_sp2dn->get_wss_internal();
    wss_persistent += utils::round_up(op_X_sp2dn->get_wss_persistent(), ator_ws_persistent->get_align());
    wss_tmp_preprocess_overlap = std::max(wss_tmp_preprocess_overlap, op_X_sp2dn->get_wss_tmp_preprocess());
    wss_tmp_perform_overlap = std::max(wss_tmp_perform_overlap, op_X_sp2dn->get_wss_tmp_perform());

    op_d_trsm.set_config(cfg.cfg_trsm);
    op_d_trsm.set_handles(q, handle_spblas, handle_dnblas);
    op_d_trsm.set_matrix_h_L(h_L_to_use);
    op_d_trsm.set_matrix_d_X(&d_X_dn);
    op_d_trsm.set_X_pattern(&h_X_sp);
    op_d_trsm.setup();
    wss_internal += op_d_trsm.get_wss_internal();
    wss_persistent += utils::round_up(op_d_trsm.get_wss_persistent(), ator_ws_persistent->get_align());
    wss_tmp_preprocess_overlap = std::max(wss_tmp_preprocess_overlap, op_d_trsm.get_wss_tmp_preprocess());
    wss_tmp_perform_overlap = std::max(wss_tmp_perform_overlap, op_d_trsm.get_wss_tmp_perform());

    op_d_herk.set_config(cfg.cfg_herk);
    op_d_herk.set_handles(q, handle_dnblas);
    op_d_herk.set_matrix_d_A(&d_X_dn);
    op_d_herk.set_matrix_d_C(&d_sc_tmp1_x);
    op_d_herk.set_h_A_pattern(&h_X_sp);
    op_d_herk.set_coefficients(-alpha, 0 * alpha); // beta=0, because I assume A22=0
    op_d_herk.set_mode(math::blas::herk_mode::AhA);
    op_d_herk.setup();
    wss_tmp_perform_overlap = std::max(wss_tmp_perform_overlap, op_d_herk.get_wss_tmp_perform());

    op_d_sc_trans = convert_ddnx_ddny<T>::make();
    op_d_sc_trans->set_handles(q, handle_dnblas);
    op_d_sc_trans->set_matrix_src(&d_sc_tmp1_x);
    op_d_sc_trans->set_matrix_dst(&d_sc_tmp2_y);
    op_d_sc_trans->setup();
    wss_tmp_perform_overlap = std::max(wss_tmp_perform_overlap, op_d_sc_trans->get_wss_tmp_perform());

    op_d_copy_sc_tmp = copy_ddnx_ddnx<T>::make();
    op_d_copy_sc_tmp->set_handles(q);
    op_d_copy_sc_tmp->set_matrix_src(&d_sc_tmp2_x);
    op_d_copy_sc_tmp->set_matrix_dst(&d_sc_tmp1_x);
    op_d_copy_sc_tmp->set_uplo(change_uplo(d_sc->prop.uplo));
    op_d_copy_sc_tmp->setup();
    wss_tmp_perform_overlap = std::max(wss_tmp_perform_overlap, op_d_copy_sc_tmp->get_wss_tmp_perform());

    op_d_perm_sc = permute_ddnx_ddnx<T,I>::make();
    op_d_perm_sc->set_handles(q);
    op_d_perm_sc->set_matrix_src(&d_sc_tmp1_x);
    op_d_perm_sc->set_matrix_dst(&d_sc_tmp2_x);
    op_d_perm_sc->set_perm_rows(&d_perm_to_sort_back_sc);
    op_d_perm_sc->set_perm_cols(&d_perm_to_sort_back_sc);
    op_d_perm_sc->setup();
    wss_tmp_perform_overlap = std::max(wss_tmp_perform_overlap, op_d_perm_sc->get_wss_tmp_perform());

    op_d_copy_sc_final = copy_ddnx_ddnx<T>::make();
    op_d_copy_sc_final->set_handles(q);
    op_d_copy_sc_final->set_matrix_src(&d_sc_tmp2_x);
    op_d_copy_sc_final->set_matrix_dst(d_sc);
    op_d_copy_sc_final->set_uplo(d_sc->prop.uplo);
    op_d_copy_sc_final->setup();
    wss_tmp_perform_overlap = std::max(wss_tmp_perform_overlap, op_d_copy_sc_final->get_wss_tmp_perform());

    wss_tmp_preprocess_linear = ((wss_tmp_preprocess_linear - 1) / ator_ws_tmp_linear->get_align() + 1) * ator_ws_tmp_linear->get_align();
    wss_tmp_perform_linear = ((wss_tmp_perform_linear - 1) / ator_ws_tmp_linear->get_align() + 1) * ator_ws_tmp_linear->get_align();

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
size_t sc_symm_hcsx_ddny_tria<T,I>::get_wss_internal()
{
    if(!called_setup) eslog::error("setup was not called\n");

    return wss_internal;
}



template<typename T, typename I>
size_t sc_symm_hcsx_ddny_tria<T,I>::get_wss_persistent()
{
    if(!called_setup) eslog::error("setup was not called\n");

    return wss_persistent;
}



template<typename T, typename I>
size_t sc_symm_hcsx_ddny_tria<T,I>::get_wss_tmp_preprocess()
{
    if(!called_setup) eslog::error("setup was not called\n");

    return wss_tmp_preprocess;
}



template<typename T, typename I>
size_t sc_symm_hcsx_ddny_tria<T,I>::get_wss_tmp_perform()
{
    if(!called_setup) eslog::error("setup was not called\n");

    return wss_tmp_perform;
}



template<typename T, typename I>
void sc_symm_hcsx_ddny_tria<T,I>::set_ws_persistent(void * ws_persistent_)
{
    if(ws_persistent_ == nullptr && wss_persistent > 0) eslog::error("persistent workspace is null\n");
    if(ws_persistent != nullptr) eslog::error("cannot re-set persistent workspace\n");

    ws_persistent = ws_persistent_;
}



template<typename T, typename I>
void sc_symm_hcsx_ddny_tria<T,I>::preprocess_submit(void * ws_tmp)
{
    stacktimer::push("sc_symm_hcsx_ddny_tria::preprocess_submit");

    if(!called_setup) eslog::error("setup has not been called\n");
    if(called_preprocess) eslog::error("preprocess has already been called\n");
    if(ws_tmp == nullptr && wss_tmp_preprocess > 0) eslog::error("temporary workspace is null\n");

    ator_ws_persistent->set(ws_persistent, wss_persistent);

    ator_ws_tmp_linear->set(ws_tmp, wss_tmp_preprocess_linear);
    ator_ws_tmp_overlap->set((char*)ws_tmp + wss_tmp_preprocess_linear, wss_tmp_preprocess_overlap);

    d_perm_to_sort_back_sc.alloc();

    gpu::mgm::copy_submit(q, h_perm_to_sort_back_sc, d_perm_to_sort_back_sc);

    op_X_sp2dn->set_ws_persistent(ator_ws_persistent->alloc(op_X_sp2dn->get_wss_persistent()));
    op_X_sp2dn->preprocess_submit(ator_ws_tmp_overlap->alloc(op_X_sp2dn->get_wss_tmp_preprocess()));

    op_d_trsm.set_ws_persistent(ator_ws_persistent->alloc(op_d_trsm.get_wss_persistent()));
    op_d_trsm.preprocess_submit(ator_ws_tmp_overlap->alloc(op_d_trsm.get_wss_tmp_preprocess()));

    ator_ws_tmp_linear->unset();
    ator_ws_tmp_overlap->unset();

    stacktimer::pop();

    called_preprocess = true;
}



template<typename T, typename I>
void sc_symm_hcsx_ddny_tria<T,I>::perform_submit(void * ws_tmp)
{
    stacktimer::push("sc_symm_hcsx_ddny_tria::perform_submit");

    if(!called_preprocess) eslog::error("preprocess has not been called\n");
    if(ws_tmp == nullptr && wss_tmp_perform > 0) eslog::error("temporary workspace is null\n");

    ator_ws_tmp_linear->set(ws_tmp, wss_tmp_perform_linear);
    ator_ws_tmp_overlap->set((char*)ws_tmp + wss_tmp_perform_linear, wss_tmp_perform_overlap);

    stacktimer::push("extract_factors");
    if(solver_factor_uplo == 'U') {
        Matrix_CSR<T,I,gpu::mgm::Ah> h_factor_old = MatrixCsxView_new<T,I>::template to_old<gpu::mgm::Ah>(h_factor_U_row);
        h_A11_solver->getFactorU(h_factor_old, false, true);
    }
    if(solver_factor_uplo == 'L') {
        Matrix_CSR<T,I,gpu::mgm::Ah> h_factor_old = MatrixCsxView_new<T,I>::template to_old<gpu::mgm::Ah>(h_factor_L_row);
        h_A11_solver->getFactorL(h_factor_old, false, true);
    }
    stacktimer::pop();

    if(need_reorder_factor_L2U) {
        op_L2U.perform_values();
    }
    if(need_reorder_factor_U2L) {
        op_U2L.perform_values();
    }

    d_X_sp.alloc();
    d_X_dn.alloc();
    d_sc_tmp1_x.alloc();
    d_sc_tmp2_x.alloc();
    d_sc_tmp1_y = d_sc_tmp1_x.get_transposed_reordered_view();
    d_sc_tmp2_y = d_sc_tmp2_x.get_transposed_reordered_view();

    gpu::mgm::copy_submit(q, h_X_sp, d_X_sp);

    op_X_sp2dn->perform_submit(ator_ws_tmp_overlap->alloc(op_X_sp2dn->get_wss_tmp_perform()));

    op_d_trsm.perform_submit(ator_ws_tmp_overlap->alloc(op_d_trsm.get_wss_tmp_perform()));

    op_d_herk.perform_submit(ator_ws_tmp_overlap->alloc(op_d_herk.get_wss_tmp_perform()));

    op_d_sc_trans->perform_submit(ator_ws_tmp_overlap->alloc(op_d_sc_trans->get_wss_tmp_perform()));

    op_d_copy_sc_tmp->perform_submit(ator_ws_tmp_overlap->alloc(op_d_copy_sc_tmp->get_wss_tmp_perform()));

    op_d_perm_sc->perform_submit(ator_ws_tmp_overlap->alloc(op_d_perm_sc->get_wss_tmp_perform()));

    op_d_copy_sc_final->perform_submit(ator_ws_tmp_overlap->alloc(op_d_copy_sc_final->get_wss_tmp_perform()));

    d_X_sp.free();
    d_X_dn.free();
    d_sc_tmp1_x.free();
    d_sc_tmp2_x.free();

    ator_ws_tmp_linear->unset();
    ator_ws_tmp_overlap->unset();

    stacktimer::pop();
}



#define INSTANTIATE_T_I(T,I) \
template class sc_symm_hcsx_ddny_tria<T,I>;

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
