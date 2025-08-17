
#include "gpu/operations/schur_hcsx_ddny.tria.h"

#include "config/ecf/operations/gpu_schur_hcsx_ddny.tria.h"
#include "esinfo/ecfinfo.h"
#include "math/wrappers/math.spsolver.h"
#include "math/operations/solver_csx.h"
#include "math/operations/permute_csx_csx.h"
#include "math/operations/pivots_trails_csx.h"
#include "math/operations/sorting_permutation.h"
#include "math/operations/quadrisect_csx_csy.h"
#include "gpu/operations/trsm_hcsx_ddny_tri.h"
#include "gpu/operations/herk_ddnx_ddny_tri.h"
#include "gpu/operations/convert_ddnx_ddny.h"
#include "gpu/operations/copy_ddnx_ddnx.h"
#include "gpu/operations/permute_ddnx_ddnx.h"
#include "esinfo/meshinfo.h"
#include "basis/utilities/stacktimer.h"



namespace espreso {
namespace gpu {
namespace operations {



template<typename T, typename I>
struct schur_hcsx_ddny_tria_data
{
    struct config
    {
        char order_X = '_';
        char order_L = '_';
        typename math::operations::solver_csx<T,I>::implementation_selector a11_solver_impl = math::operations::solver_csx<T,I>::implementation_selector::autoselect;
    };
    config cfg;
    MatrixCsxData_new<T,I> h_U_data;
    MatrixCsxData_new<T,I> h_L_data;
    MatrixCsxView_new<T,I> h_L_row;
    MatrixCsxView_new<T,I> h_L_col;
    MatrixCsxView_new<T,I> h_U_row;
    MatrixCsxView_new<T,I> h_U_col;
    MatrixCsxView_new<T,I> * h_L_to_use = nullptr;
    PermutationData_new<I> h_perm_to_sort_A12_cols;
    PermutationView_new<I> h_perm_to_sort_back_sc;
    PermutationData_new<I> d_perm_to_sort_back_sc;
    PermutationData_new<I> h_perm_fillreduce;
    MatrixCsxData_new<T,I> h_X_sp;
    MatrixCsxData_new<T,I> d_X_sp;
    MatrixDenseData_new<T> d_X_dn;
    MatrixDenseData_new<T> d_sc_tmp1_x; // same order as sc
    MatrixDenseData_new<T> d_sc_tmp2_x; // same order as sc
    MatrixDenseView_new<T> d_sc_tmp1_y; // different order as sc
    MatrixDenseView_new<T> d_sc_tmp2_y; // different order as sc
    std::unique_ptr<math::operations::solver_csx<T,I>> op_h_A11_solver;
    std::unique_ptr<convert_dcsx_ddny<T,I>> op_X_sp2dn;
    math::operations::convert_csx_csy_map<T,I> op_L2U;
    math::operations::convert_csx_csy_map<T,I> op_U2L;
    trsm_hcsx_ddny_tri<T,I> op_d_trsm;
    herk_ddnx_ddny_tri<T,I> op_d_herk;
    std::unique_ptr<convert_ddnx_ddny<T>> op_d_sc_trans;
    std::unique_ptr<copy_ddnx_ddnx<T>> op_d_copy_sc_tmp;
    std::unique_ptr<permute_ddnx_ddnx<T,I>> op_d_perm_sc;
    std::unique_ptr<copy_ddnx_ddnx<T>> op_d_copy_sc_final;
    math::operations::quadrisect_csx_csy<T,I> op_split;
    MatrixCsxData_new<T,I> sub_h_A11;
    MatrixCsxData_new<T,I> sub_h_A12;
    MatrixCsxData_new<T,I> sub_h_A21;
    MatrixCsxData_new<T,I> sub_h_A22;
    char solver_factor_uplo = '_';
    bool need_reorder_factor_L2U = false;
    bool need_reorder_factor_U2L = false;
};



template<typename T, typename I>
void setup_config(typename schur_hcsx_ddny_tria_data<T,I>::config & cfg)
{
    using solver_impl_t = typename math::operations::solver_csx<T,I>::implementation_selector;
    using ecf_config = GpuSchurHcsxDdnyTriaConfig;
    const ecf_config & ecf = info::ecf->operations.gpu_schur_hcsx_ddny_tria;

    switch(ecf.order_X) {
        case ecf_config::MATRIX_ORDER::AUTO:      cfg.order_X = 'R'; break;
        case ecf_config::MATRIX_ORDER::ROW_MAJOR: cfg.order_X = 'R'; break;
        case ecf_config::MATRIX_ORDER::COL_MAJOR: cfg.order_X = 'C'; break;
    }

    switch(ecf.order_L) {
        case ecf_config::MATRIX_ORDER::AUTO:      cfg.order_L = 'C'; break;
        case ecf_config::MATRIX_ORDER::ROW_MAJOR: cfg.order_L = 'R'; break;
        case ecf_config::MATRIX_ORDER::COL_MAJOR: cfg.order_L = 'C'; break;
    }

    switch(ecf.sparse_solver_impl) {
        case ecf_config::SPARSE_SOLVER_IMPL::AUTO:         cfg.a11_solver_impl = solver_impl_t::autoselect;   break;
        case ecf_config::SPARSE_SOLVER_IMPL::MKLPARDISO:   cfg.a11_solver_impl = solver_impl_t::mklpardiso;   break;
        case ecf_config::SPARSE_SOLVER_IMPL::SUITESPARSE:  cfg.a11_solver_impl = solver_impl_t::suitesparse;  break;
        case ecf_config::SPARSE_SOLVER_IMPL::MUMPS:        cfg.a11_solver_impl = solver_impl_t::mumps;        break;
        case ecf_config::SPARSE_SOLVER_IMPL::STRUMPACK:    cfg.a11_solver_impl = solver_impl_t::strumpack;    break;
        case ecf_config::SPARSE_SOLVER_IMPL::PASTIX:       cfg.a11_solver_impl = solver_impl_t::pastix;       break;
        case ecf_config::SPARSE_SOLVER_IMPL::SUPERLU_DIST: cfg.a11_solver_impl = solver_impl_t::superlu_dist; break;
    }
}



template<typename T, typename I>
schur_hcsx_ddny_tria<T,I>::schur_hcsx_ddny_tria() = default;



template<typename T, typename I>
schur_hcsx_ddny_tria<T,I>::~schur_hcsx_ddny_tria() = default;



template<typename T, typename I>
void schur_hcsx_ddny_tria<T,I>::internal_setup()
{
    if(!is_matrix_hermitian) eslog::error("for now only hermitian matrices are supported\n");
    if(h_A12 == nullptr) eslog::error("the upper part of matrix must be stored\n");
    if(h_A11->prop.uplo != 'U') eslog::error("the upper part of A11 must be stored\n");
    if(h_A22 != nullptr) eslog::error("for now I assume A22=O\n");

    data = std::make_unique<schur_hcsx_ddny_tria_data<T,I>>();
    setup_config<T,I>(data->cfg);

    data->solver_factor_uplo = h_A11->prop.uplo;

    if(called_set_matrix == '1') {
        data->op_split.set_matrix_src(h_A);
        data->op_split.set_matrices_dst(&data->sub_h_A11, &data->sub_h_A12, &data->sub_h_A21, &data->sub_h_A22);
        data->op_split.set_bounds(size_A11, size_A11);
        data->op_split.setup();

        data->sub_h_A11.set(size_A11, size_A11, data->op_split.get_output_matrix_11_nnz(), h_A->order, AllocatorHostPinned_new::get_singleton());
        data->sub_h_A12.set(size_A11, size_sc,  data->op_split.get_output_matrix_12_nnz(), h_A->order, AllocatorHostPinned_new::get_singleton());
        data->sub_h_A21.set(size_sc,  size_A11, data->op_split.get_output_matrix_21_nnz(), h_A->order, AllocatorHostPinned_new::get_singleton());
        data->sub_h_A22.set(size_sc,  size_sc,  data->op_split.get_output_matrix_22_nnz(), h_A->order, AllocatorHostPinned_new::get_singleton());

        data->sub_h_A11.prop = h_A->prop;
        data->sub_h_A22.prop = h_A->prop;

        data->sub_h_A11.alloc();
        data->sub_h_A12.alloc();
        data->sub_h_A21.alloc();
        data->sub_h_A22.alloc();

        data->op_split.perform();

        h_A11 = &data->sub_h_A11;
        if(h_A->prop.uplo != 'L') h_A12 = &data->sub_h_A12;
        if(h_A->prop.uplo != 'U') h_A21 = &data->sub_h_A21;
        h_A22 = &data->sub_h_A22;
    }

    data->op_h_A11_solver = math::operations::solver_csx<T,I>::make(data->cfg.a11_solver_impl, &h_A11->prop, true, need_solve_A11);
    data->op_h_A11_solver->set_matrix_A(h_A11);
    data->op_h_A11_solver->set_needs(true, need_solve_A11);
    data->op_h_A11_solver->factorize_symbolic();

    size_t factor_nnz = data->op_h_A11_solver->get_factor_nnz();

    data->h_U_data.set(size_A11, size_A11, factor_nnz, h_A11->order, AllocatorHostPinned_new::get_singleton());
    data->h_U_data.prop.uplo = 'U';
    data->h_U_data.prop.diag = 'N';

    data->h_L_data.set(size_A11, size_A11, factor_nnz, h_A11->order, AllocatorHostPinned_new::get_singleton());
    data->h_L_data.prop.uplo = 'L';
    data->h_L_data.prop.diag = 'N';

    data->need_reorder_factor_L2U = ((data->solver_factor_uplo == 'L') && (data->cfg.order_L == 'C'));
    data->need_reorder_factor_U2L = ((data->solver_factor_uplo == 'U') && (data->cfg.order_L == 'R'));

    if(data->solver_factor_uplo == 'U' || data->need_reorder_factor_L2U) {
        data->h_U_data.alloc();
    }
    if(data->solver_factor_uplo == 'L' || data->need_reorder_factor_U2L) {
        data->h_L_data.alloc();
    }

    if(h_A11->prop.uplo == 'L') {
        data->op_h_A11_solver->get_factor_L(data->h_L_data, true, false);
    }
    if(h_A11->prop.uplo == 'U') {
        data->op_h_A11_solver->get_factor_U(data->h_U_data, true, false);
    }

    if(h_A11->order == 'R') {
        data->h_L_row = data->h_L_data;
        data->h_U_row = data->h_U_data;
        data->h_L_col = data->h_U_row.get_transposed_reordered_view();
        data->h_U_col = data->h_L_row.get_transposed_reordered_view();
    }
    if(h_A11->order == 'C') {
        data->h_L_col = data->h_L_data;
        data->h_U_col = data->h_U_data;
        data->h_L_row = data->h_U_col.get_transposed_reordered_view();
        data->h_U_row = data->h_L_col.get_transposed_reordered_view();
    }

    if(data->need_reorder_factor_L2U) {
        data->op_L2U.set_matrix_src(&data->h_L_row);
        data->op_L2U.set_matrix_dst(&data->h_L_col);
        data->op_L2U.perform_pattern();
    }
    if(data->need_reorder_factor_U2L) {
        data->op_U2L.set_matrix_src(&data->h_U_row);
        data->op_U2L.set_matrix_dst(&data->h_U_col);
        data->op_U2L.perform_pattern();
    }

    if(data->cfg.order_L == 'R') data->h_L_to_use = &data->h_L_row;
    if(data->cfg.order_L == 'C') data->h_L_to_use = &data->h_L_col;

    data->h_X_sp.set(h_A12->nrows, h_A12->ncols, h_A12->nnz, h_A12->order, AllocatorHostPinned_new::get_singleton());
    data->h_X_sp.alloc();

    {
        MatrixCsxData_new<T,I> h_X_sp_tmp;
        h_X_sp_tmp.set(data->h_X_sp.nrows, data->h_X_sp.ncols, data->h_X_sp.nnz, data->h_X_sp.order, AllocatorCPU_new::get_singleton());
        h_X_sp_tmp.alloc();

        {
            data->h_perm_fillreduce.set(h_A12->nrows, AllocatorCPU_new::get_singleton());
            data->h_perm_fillreduce.alloc();
            data->op_h_A11_solver->get_permutation(data->h_perm_fillreduce);

            math::operations::permute_csx_csx<T,I>::do_all(h_A12, &h_X_sp_tmp, &data->h_perm_fillreduce, nullptr);
        }

        {
            VectorDenseData_new<I> colpivots;
            colpivots.set(h_A12->ncols, AllocatorCPU_new::get_singleton());
            colpivots.alloc();

            math::operations::pivots_trails_csx<T,I>::do_all(&h_X_sp_tmp, &colpivots, 'C', 'P', '_');

            data->h_perm_to_sort_A12_cols.set(h_A12->ncols, AllocatorHostPinned_new::get_singleton());
            data->h_perm_to_sort_A12_cols.alloc();

            data->h_perm_to_sort_back_sc = data->h_perm_to_sort_A12_cols.get_inverse_view();

            data->d_perm_to_sort_back_sc.set(data->h_perm_to_sort_back_sc.size, ator_ws_persistent.get());
            wss_persistent += utils::round_up(data->d_perm_to_sort_back_sc.get_memory_impact(), ator_ws_persistent->get_align());

            math::operations::sorting_permutation<I,I>::do_all(&colpivots, &data->h_perm_to_sort_A12_cols);
        }

        math::operations::permute_csx_csx<T,I>::do_all(&h_X_sp_tmp, &data->h_X_sp, nullptr, &data->h_perm_to_sort_A12_cols);
    }

    data->d_X_sp.set(data->h_X_sp.nrows, data->h_X_sp.ncols, data->h_X_sp.nnz, data->h_X_sp.order, ator_ws_tmp_linear.get());
    wss_tmp_perform_linear += data->d_X_sp.get_memory_impact();

    data->d_X_dn.set(data->h_X_sp.nrows, data->h_X_sp.ncols, data->cfg.order_X, ator_ws_tmp_linear.get());
    wss_tmp_perform_linear += data->d_X_dn.get_memory_impact();

    data->d_sc_tmp1_x.set(d_sc->nrows, d_sc->ncols, d_sc->order, ator_ws_tmp_linear.get());
    data->d_sc_tmp1_x.prop.uplo = d_sc->prop.uplo;
    wss_tmp_perform_linear += data->d_sc_tmp1_x.get_memory_impact();
    data->d_sc_tmp1_y = data->d_sc_tmp1_x.get_transposed_reordered_view();

    data->d_sc_tmp2_x.set(d_sc->nrows, d_sc->ncols, d_sc->order, ator_ws_tmp_linear.get());
    data->d_sc_tmp2_x.prop.uplo = d_sc->prop.uplo;
    wss_tmp_perform_linear += data->d_sc_tmp2_x.get_memory_impact();
    data->d_sc_tmp2_y = data->d_sc_tmp2_x.get_transposed_reordered_view();

    data->op_X_sp2dn = convert_dcsx_ddny<T,I>::make();
    data->op_X_sp2dn->set_handles(q, handle_spblas);
    data->op_X_sp2dn->set_matrix_src(&data->d_X_sp);
    data->op_X_sp2dn->set_matrix_dst(&data->d_X_dn);
    data->op_X_sp2dn->setup();
    wss_internal += data->op_X_sp2dn->get_wss_internal();
    wss_persistent += utils::round_up(data->op_X_sp2dn->get_wss_persistent(), ator_ws_persistent->get_align());
    wss_tmp_preprocess_overlap = std::max(wss_tmp_preprocess_overlap, data->op_X_sp2dn->get_wss_tmp_preprocess());
    wss_tmp_perform_overlap = std::max(wss_tmp_perform_overlap, data->op_X_sp2dn->get_wss_tmp_perform());

    data->op_d_trsm.set_handles(q, handle_spblas, handle_dnblas);
    data->op_d_trsm.set_matrix_h_L(data->h_L_to_use);
    data->op_d_trsm.set_matrix_d_X(&data->d_X_dn);
    data->op_d_trsm.set_X_pattern(&data->h_X_sp);
    data->op_d_trsm.setup();
    wss_internal += data->op_d_trsm.get_wss_internal();
    wss_persistent += utils::round_up(data->op_d_trsm.get_wss_persistent(), ator_ws_persistent->get_align());
    wss_tmp_preprocess_overlap = std::max(wss_tmp_preprocess_overlap, data->op_d_trsm.get_wss_tmp_preprocess());
    wss_tmp_perform_overlap = std::max(wss_tmp_perform_overlap, data->op_d_trsm.get_wss_tmp_perform());

    data->op_d_herk.set_handles(q, handle_dnblas);
    data->op_d_herk.set_matrix_d_A(&data->d_X_dn);
    data->op_d_herk.set_matrix_d_C(&data->d_sc_tmp1_x);
    data->op_d_herk.set_h_A_pattern(&data->h_X_sp);
    data->op_d_herk.set_coefficients(-alpha, 0 * alpha); // beta=0, because I assume A22=0
    data->op_d_herk.set_mode(math::blas::herk_mode::AhA);
    data->op_d_herk.setup();
    wss_tmp_perform_overlap = std::max(wss_tmp_perform_overlap, data->op_d_herk.get_wss_tmp_perform());

    data->op_d_sc_trans = convert_ddnx_ddny<T>::make();
    data->op_d_sc_trans->set_handles(q, handle_dnblas);
    data->op_d_sc_trans->set_matrix_src(&data->d_sc_tmp1_x);
    data->op_d_sc_trans->set_matrix_dst(&data->d_sc_tmp2_y);
    data->op_d_sc_trans->setup();
    wss_tmp_perform_overlap = std::max(wss_tmp_perform_overlap, data->op_d_sc_trans->get_wss_tmp_perform());

    data->op_d_copy_sc_tmp = copy_ddnx_ddnx<T>::make();
    data->op_d_copy_sc_tmp->set_handles(q);
    data->op_d_copy_sc_tmp->set_matrix_src(&data->d_sc_tmp2_x);
    data->op_d_copy_sc_tmp->set_matrix_dst(&data->d_sc_tmp1_x);
    data->op_d_copy_sc_tmp->set_uplo(change_uplo(d_sc->prop.uplo));
    data->op_d_copy_sc_tmp->setup();
    wss_tmp_perform_overlap = std::max(wss_tmp_perform_overlap, data->op_d_copy_sc_tmp->get_wss_tmp_perform());

    data->op_d_perm_sc = permute_ddnx_ddnx<T,I>::make();
    data->op_d_perm_sc->set_handles(q);
    data->op_d_perm_sc->set_matrix_src(&data->d_sc_tmp1_x);
    data->op_d_perm_sc->set_matrix_dst(&data->d_sc_tmp2_x);
    data->op_d_perm_sc->set_perm_rows(&data->d_perm_to_sort_back_sc);
    data->op_d_perm_sc->set_perm_cols(&data->d_perm_to_sort_back_sc);
    data->op_d_perm_sc->setup();
    wss_tmp_perform_overlap = std::max(wss_tmp_perform_overlap, data->op_d_perm_sc->get_wss_tmp_perform());

    data->op_d_copy_sc_final = copy_ddnx_ddnx<T>::make();
    data->op_d_copy_sc_final->set_handles(q);
    data->op_d_copy_sc_final->set_matrix_src(&data->d_sc_tmp2_x);
    data->op_d_copy_sc_final->set_matrix_dst(d_sc);
    data->op_d_copy_sc_final->set_uplo(d_sc->prop.uplo);
    data->op_d_copy_sc_final->setup();
    wss_tmp_perform_overlap = std::max(wss_tmp_perform_overlap, data->op_d_copy_sc_final->get_wss_tmp_perform());

    wss_tmp_preprocess_linear = ((wss_tmp_preprocess_linear - 1) / ator_ws_tmp_linear->get_align() + 1) * ator_ws_tmp_linear->get_align();
    wss_tmp_perform_linear = ((wss_tmp_perform_linear - 1) / ator_ws_tmp_linear->get_align() + 1) * ator_ws_tmp_linear->get_align();

    wss_tmp_preprocess = wss_tmp_preprocess_linear + wss_tmp_preprocess_overlap;
    wss_tmp_perform = wss_tmp_perform_linear + wss_tmp_perform_overlap;

    stacktimer::info("wss_internal       %zu", wss_internal);
    stacktimer::info("wss_persistent     %zu", wss_persistent);
    stacktimer::info("wss_tmp_preprocess %zu", wss_tmp_preprocess);
    stacktimer::info("wss_tmp_perform    %zu", wss_tmp_perform);

    called_setup = true;
}



template<typename T, typename I>
void schur_hcsx_ddny_tria<T,I>::internal_preprocess_submit()
{
    data->d_perm_to_sort_back_sc.alloc();

    gpu::mgm::copy_submit(q, data->h_perm_to_sort_back_sc, data->d_perm_to_sort_back_sc);

    data->op_X_sp2dn->set_ws_persistent(ator_ws_persistent->alloc(data->op_X_sp2dn->get_wss_persistent()));
    data->op_X_sp2dn->preprocess_submit(ator_ws_tmp_overlap->alloc(data->op_X_sp2dn->get_wss_tmp_preprocess()));

    data->op_d_trsm.set_ws_persistent(ator_ws_persistent->alloc(data->op_d_trsm.get_wss_persistent()));
    data->op_d_trsm.preprocess_submit(ator_ws_tmp_overlap->alloc(data->op_d_trsm.get_wss_tmp_preprocess()));
}



template<typename T, typename I>
void schur_hcsx_ddny_tria<T,I>::internal_perform_1_submit()
{
    if(called_set_matrix == '1') {
        data->op_split.perform();
    }

    stacktimer::push("numerical_factorization");
    data->op_h_A11_solver->factorize_numeric();
    stacktimer::pop();
}



template<typename T, typename I>
void schur_hcsx_ddny_tria<T,I>::internal_perform_2_submit()
{
    stacktimer::push("extract_factors");
    if(data->solver_factor_uplo == 'U') {
        data->op_h_A11_solver->get_factor_U(data->h_U_data, false, true);
    }
    if(data->solver_factor_uplo == 'L') {
        data->op_h_A11_solver->get_factor_L(data->h_L_data, false, true);
    }
    stacktimer::pop();

    if(data->need_reorder_factor_L2U) {
        data->op_L2U.perform_values();
    }
    if(data->need_reorder_factor_U2L) {
        data->op_U2L.perform_values();
    }

    data->d_X_sp.alloc();
    data->d_X_dn.alloc();
    data->d_sc_tmp1_x.alloc();
    data->d_sc_tmp2_x.alloc();
    data->d_sc_tmp1_y = data->d_sc_tmp1_x.get_transposed_reordered_view();
    data->d_sc_tmp2_y = data->d_sc_tmp2_x.get_transposed_reordered_view();

    gpu::mgm::copy_submit(q, data->h_X_sp, data->d_X_sp);

    data->op_X_sp2dn->perform_submit(ator_ws_tmp_overlap->alloc(data->op_X_sp2dn->get_wss_tmp_perform()));

    data->op_d_trsm.perform_submit(ator_ws_tmp_overlap->alloc(data->op_d_trsm.get_wss_tmp_perform()));

    data->op_d_herk.perform_submit(ator_ws_tmp_overlap->alloc(data->op_d_herk.get_wss_tmp_perform()));

    data->op_d_sc_trans->perform_submit(ator_ws_tmp_overlap->alloc(data->op_d_sc_trans->get_wss_tmp_perform()));

    data->op_d_copy_sc_tmp->perform_submit(ator_ws_tmp_overlap->alloc(data->op_d_copy_sc_tmp->get_wss_tmp_perform()));

    data->op_d_perm_sc->perform_submit(ator_ws_tmp_overlap->alloc(data->op_d_perm_sc->get_wss_tmp_perform()));

    data->op_d_copy_sc_final->perform_submit(ator_ws_tmp_overlap->alloc(data->op_d_copy_sc_final->get_wss_tmp_perform()));

    data->d_X_sp.free();
    data->d_X_dn.free();
    data->d_sc_tmp1_x.free();
    data->d_sc_tmp2_x.free();
}



template<typename T, typename I>
void schur_hcsx_ddny_tria<T,I>::internal_solve_A11(VectorDenseView_new<T> & rhs, VectorDenseView_new<T> & sol)
{
    data->op_h_A11_solver->solve(rhs, sol);
}



#define INSTANTIATE_T_I(T,I) \
template class schur_hcsx_ddny_tria<T,I>;

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
