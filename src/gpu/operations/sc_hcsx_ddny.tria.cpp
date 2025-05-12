
#include "gpu/operations/sc_hcsx_ddny.tria.h"

#include "math/wrappers/math.spsolver.h"
#include "math/primitives_new/permutation_data_new.h"
#include "math/operations/permute_csx_csx.h"
#include "math/operations/pivots_trails_csx.h"
#include "math/operations/sorting_permutation.h"
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
struct sc_hcsx_ddny_tria_data
{
    struct config
    {
        typename trsm_hcsx_ddny_tri<T,I>::config cfg_trsm;
        typename herk_ddnx_ddny_tri<T,I>::config cfg_herk;
        char order_X = '_';
        char order_L = '_';
    };
    config cfg;
    DirectSparseSolver<T,I> h_A11_solver;
    Matrix_CSR<T,I> h_A11_old;
    MatrixCsxData_new<T,I> h_factor_U_row;
    MatrixCsxData_new<T,I> h_factor_L_row;
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
    std::unique_ptr<convert_dcsx_ddny<T,I>> op_X_sp2dn;
    math::operations::convert_csx_csy_map<T,I> op_L2U;
    math::operations::convert_csx_csy_map<T,I> op_U2L;
    trsm_hcsx_ddny_tri<T,I> op_d_trsm;
    herk_ddnx_ddny_tri<T,I> op_d_herk;
    std::unique_ptr<convert_ddnx_ddny<T>> op_d_sc_trans;
    std::unique_ptr<copy_ddnx_ddnx<T>> op_d_copy_sc_tmp;
    std::unique_ptr<permute_ddnx_ddnx<T,I>> op_d_perm_sc;
    std::unique_ptr<copy_ddnx_ddnx<T>> op_d_copy_sc_final;
    char solver_factor_uplo = '_';
    bool need_reorder_factor_L2U = false;
    bool need_reorder_factor_U2L = false;
};



template<typename T>
static void set_by_env(T & var, const char * env_var)
{
    const char * str = std::getenv(env_var);
    if(str != nullptr) {
        if constexpr(std::is_same_v<T,char>) var = *str;
        if constexpr(std::is_same_v<T,int>) var = atoi(str);
        if constexpr(std::is_same_v<T,double>) var = atof(str);
        if constexpr(std::is_same_v<T,bool>) {
            var = (strcmp(str,"1") == 0
                || strcmp(str,"Y") == 0
                || strcmp(str,"Yes") == 0
                || strcmp(str,"yes") == 0
                || strcmp(str,"T") == 0
                || strcmp(str,"True") == 0
                || strcmp(str,"true") == 0
                || strcmp(str,"On") == 0
                || strcmp(str,"on") == 0);
        }
    }
}

static void replace_if_default(char & param, char deflt)
{
    if(param == '_') {
        param = deflt;
    }
}

static void replace_if_zero(int & param, int deflt)
{
    if(param == 0) {
        param = deflt;
    }
}

template<typename T, typename I>
static void populate_config_from_env(typename sc_hcsx_ddny_tria_data<T,I>::config & cfg)
{
    cfg.cfg_trsm.partition.parameter = 0;
    cfg.cfg_herk.partition_parameter = 0;

    set_by_env(cfg.order_X,                                   "ESPRESO_CONFIG_OP_SC_HCSX_DDNY_TRIA_order_X");
    set_by_env(cfg.order_L,                                   "ESPRESO_CONFIG_OP_SC_HCSX_DDNY_TRIA_order_L");
    set_by_env(cfg.cfg_trsm.strategy,                         "ESPRESO_CONFIG_OP_SC_HCSX_DDNY_TRIA_trsm_strategy");
    set_by_env(cfg.cfg_trsm.partition.algorithm,              "ESPRESO_CONFIG_OP_SC_HCSX_DDNY_TRIA_trsm_partition_algorithm");
    set_by_env(cfg.cfg_trsm.partition.parameter,              "ESPRESO_CONFIG_OP_SC_HCSX_DDNY_TRIA_trsm_partition_parameter");
    set_by_env(cfg.cfg_trsm.splitrhs.factor_order_sp,         "ESPRESO_CONFIG_OP_SC_HCSX_DDNY_TRIA_trsm_splitrhs_factor_order_sp");
    set_by_env(cfg.cfg_trsm.splitrhs.factor_order_dn,         "ESPRESO_CONFIG_OP_SC_HCSX_DDNY_TRIA_trsm_splitrhs_factor_order_dn");
    set_by_env(cfg.cfg_trsm.splitrhs.spdn_criteria,           "ESPRESO_CONFIG_OP_SC_HCSX_DDNY_TRIA_trsm_splitrhs_spdn_criteria");
    set_by_env(cfg.cfg_trsm.splitrhs.spdn_param,              "ESPRESO_CONFIG_OP_SC_HCSX_DDNY_TRIA_trsm_splitrhs_spdn_param");
    set_by_env(cfg.cfg_trsm.splitfactor.trsm_factor_spdn,     "ESPRESO_CONFIG_OP_SC_HCSX_DDNY_TRIA_trsm_splitfactor_trsm_factor_spdn");
    set_by_env(cfg.cfg_trsm.splitfactor.trsm_factor_order,    "ESPRESO_CONFIG_OP_SC_HCSX_DDNY_TRIA_trsm_splitfactor_trsm_factor_order");
    set_by_env(cfg.cfg_trsm.splitfactor.gemm_factor_order_sp, "ESPRESO_CONFIG_OP_SC_HCSX_DDNY_TRIA_trsm_splitfactor_gemm_factor_order_sp");
    set_by_env(cfg.cfg_trsm.splitfactor.gemm_factor_order_dn, "ESPRESO_CONFIG_OP_SC_HCSX_DDNY_TRIA_trsm_splitfactor_gemm_factor_order_dn");
    set_by_env(cfg.cfg_trsm.splitfactor.gemm_factor_prune,    "ESPRESO_CONFIG_OP_SC_HCSX_DDNY_TRIA_trsm_splitfactor_gemm_factor_prune");
    set_by_env(cfg.cfg_trsm.splitfactor.gemm_spdn_criteria,   "ESPRESO_CONFIG_OP_SC_HCSX_DDNY_TRIA_trsm_splitfactor_gemm_spdn_criteria");
    set_by_env(cfg.cfg_trsm.splitfactor.gemm_spdn_param,      "ESPRESO_CONFIG_OP_SC_HCSX_DDNY_TRIA_trsm_splitfactor_gemm_spdn_param");
    set_by_env(cfg.cfg_herk.strategy,                         "ESPRESO_CONFIG_OP_SC_HCSX_DDNY_TRIA_herk_strategy");
    set_by_env(cfg.cfg_herk.partition_algorithm,              "ESPRESO_CONFIG_OP_SC_HCSX_DDNY_TRIA_herk_partition_algorithm");
    set_by_env(cfg.cfg_herk.partition_parameter,              "ESPRESO_CONFIG_OP_SC_HCSX_DDNY_TRIA_herk_partition_parameter");
}

template<typename T, typename I>
static void populate_config_replace_defaults(typename sc_hcsx_ddny_tria_data<T,I>::config & cfg)
{
    replace_if_default(cfg.cfg_trsm.strategy, 'F');
    replace_if_default(cfg.cfg_herk.strategy, 'Q');

    if(info::mesh->dimension == 2 && cfg.cfg_herk.strategy == 'Q') replace_if_zero(cfg.cfg_herk.partition_parameter, -2000);
    if(info::mesh->dimension == 2 && cfg.cfg_herk.strategy == 'T') replace_if_zero(cfg.cfg_herk.partition_parameter, -200);
    if(info::mesh->dimension == 3 && cfg.cfg_herk.strategy == 'Q') replace_if_zero(cfg.cfg_herk.partition_parameter, -1000);
    if(info::mesh->dimension == 3 && cfg.cfg_herk.strategy == 'T') replace_if_zero(cfg.cfg_herk.partition_parameter, -1000);
    replace_if_default(cfg.cfg_herk.partition_algorithm, 'U');

    if(info::mesh->dimension == 2 && cfg.cfg_trsm.strategy == 'F') replace_if_zero(cfg.cfg_trsm.partition.parameter, -1000);
    if(info::mesh->dimension == 2 && cfg.cfg_trsm.strategy == 'R') replace_if_zero(cfg.cfg_trsm.partition.parameter, 1);
    if(info::mesh->dimension == 3 && cfg.cfg_trsm.strategy == 'F') replace_if_zero(cfg.cfg_trsm.partition.parameter, -500);
    if(info::mesh->dimension == 3 && cfg.cfg_trsm.strategy == 'R') replace_if_zero(cfg.cfg_trsm.partition.parameter, -1000);
    replace_if_default(cfg.cfg_trsm.partition.algorithm, 'U');

    replace_if_default(cfg.cfg_trsm.splitrhs.spdn_criteria, 'S');
    replace_if_default(cfg.cfg_trsm.splitfactor.trsm_factor_spdn, (info::mesh->dimension == 2) ? 'S' : 'D');
    replace_if_default(cfg.cfg_trsm.splitfactor.gemm_factor_prune, 'R');
    replace_if_default(cfg.cfg_trsm.splitfactor.gemm_spdn_criteria, (info::mesh->dimension == 3 && cfg.cfg_trsm.splitfactor.gemm_factor_prune == 'R') ? 'D' : 'S');

    replace_if_default(cfg.order_X, 'R');
    replace_if_default(cfg.cfg_trsm.splitrhs.factor_order_sp, 'R');
    replace_if_default(cfg.cfg_trsm.splitrhs.factor_order_dn, 'R');
    replace_if_default(cfg.cfg_trsm.splitfactor.trsm_factor_order, 'R');
    replace_if_default(cfg.cfg_trsm.splitfactor.gemm_factor_order_sp, 'R');
    replace_if_default(cfg.cfg_trsm.splitfactor.gemm_factor_order_dn, 'R');

    replace_if_default(cfg.order_L, 'C');
}

template<typename T, typename I>
static void populate_config(typename sc_hcsx_ddny_tria_data<T,I>::config & cfg)
{
    populate_config_from_env<T,I>(cfg);

    populate_config_replace_defaults<T,I>(cfg);
}



template<typename T, typename I>
sc_hcsx_ddny_tria<T,I>::sc_hcsx_ddny_tria() = default;



template<typename T, typename I>
sc_hcsx_ddny_tria<T,I>::~sc_hcsx_ddny_tria() = default;



template<typename T, typename I>
void sc_hcsx_ddny_tria<T,I>::internal_setup()
{
    if(called_set_matrix != '4') eslog::error("support only 4 small matrices on input for now\n");
    if(!is_matrix_hermitian) eslog::error("for now only hermitian matrices are supported\n");
    if(h_A12 == nullptr) eslog::error("the upper part of matrix must be stored\n");
    if(h_A11->prop.uplo != 'U') eslog::error("the upper part of A11 must be stored\n");
    if(h_A22 != nullptr) eslog::error("for now I assume A22=O\n");

    data = std::make_unique<sc_hcsx_ddny_tria_data<T,I>>();

    populate_config<T,I>(data->cfg);

    data->solver_factor_uplo = '_';
    if(DirectSparseSolver<T,I>::factorsSymmetry() == Solver_Factors::HERMITIAN_LOWER) data->solver_factor_uplo = 'L';
    if(DirectSparseSolver<T,I>::factorsSymmetry() == Solver_Factors::HERMITIAN_UPPER) data->solver_factor_uplo = 'U';
    if(data->solver_factor_uplo == '_') eslog::error("wrong sparse solver, must be symmetric\n");
    if(!DirectSparseSolver<T,I>::provideFactors()) eslog::error("wrong sparse solver, must provide factors\n");

    stacktimer::push("symbolic_factorization");
    data->h_A11_old = MatrixCsxView_new<T,I>::template to_old<cpu_allocator>(*h_A11);
    data->h_A11_solver.commit(data->h_A11_old);
    data->h_A11_solver.symbolicFactorization();
    stacktimer::pop();

    I factor_nnz = data->h_A11_solver.getFactorNnz();

    data->h_factor_U_row.set(size_A11, size_A11, factor_nnz, 'R', AllocatorHostPinned_new::get_singleton());
    data->h_factor_U_row.prop.uplo = 'U';
    data->h_factor_U_row.prop.diag = 'N';

    data->h_factor_L_row.set(size_A11, size_A11, factor_nnz, 'R', AllocatorHostPinned_new::get_singleton());
    data->h_factor_L_row.prop.uplo = 'L';
    data->h_factor_L_row.prop.diag = 'N';

    data->need_reorder_factor_L2U = ((data->solver_factor_uplo == 'L') && (data->cfg.order_L == 'C'));
    data->need_reorder_factor_U2L = ((data->solver_factor_uplo == 'U') && (data->cfg.order_L == 'R'));

    stacktimer::push("alloc_host_factors");
    if(data->solver_factor_uplo == 'U' || data->need_reorder_factor_L2U) {
        data->h_factor_U_row.alloc();
    }
    if(data->solver_factor_uplo == 'L' || data->need_reorder_factor_U2L) {
        data->h_factor_L_row.alloc();
    }
    stacktimer::pop();

    stacktimer::push("extract_factors");
    if(data->solver_factor_uplo == 'U') {
        Matrix_CSR<T,I,gpu::mgm::Ah> h_factor_old_U = MatrixCsxView_new<T,I>::template to_old<gpu::mgm::Ah>(data->h_factor_U_row);
        data->h_A11_solver.getFactorU(h_factor_old_U, true, false);
    }
    if(data->solver_factor_uplo == 'L') {
        Matrix_CSR<T,I,gpu::mgm::Ah> h_factor_old_L = MatrixCsxView_new<T,I>::template to_old<gpu::mgm::Ah>(data->h_factor_L_row);
        data->h_A11_solver.getFactorL(h_factor_old_L, true, false);
    }
    stacktimer::pop();

    data->h_L_row = data->h_factor_L_row;
    data->h_L_col = data->h_factor_U_row.get_transposed_reordered_view();
    data->h_U_row = data->h_factor_U_row;
    data->h_U_col = data->h_factor_L_row.get_transposed_reordered_view();

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
            Permutation<I> h_perm_fillreduce_old = PermutationView_new<I>::template to_old<cpu_allocator>(data->h_perm_fillreduce);
            data->h_A11_solver.getPermutation(h_perm_fillreduce_old);

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

    data->op_d_trsm.set_config(data->cfg.cfg_trsm);
    data->op_d_trsm.set_handles(q, handle_spblas, handle_dnblas);
    data->op_d_trsm.set_matrix_h_L(data->h_L_to_use);
    data->op_d_trsm.set_matrix_d_X(&data->d_X_dn);
    data->op_d_trsm.set_X_pattern(&data->h_X_sp);
    data->op_d_trsm.setup();
    wss_internal += data->op_d_trsm.get_wss_internal();
    wss_persistent += utils::round_up(data->op_d_trsm.get_wss_persistent(), ator_ws_persistent->get_align());
    wss_tmp_preprocess_overlap = std::max(wss_tmp_preprocess_overlap, data->op_d_trsm.get_wss_tmp_preprocess());
    wss_tmp_perform_overlap = std::max(wss_tmp_perform_overlap, data->op_d_trsm.get_wss_tmp_perform());

    data->op_d_herk.set_config(data->cfg.cfg_herk);
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
void sc_hcsx_ddny_tria<T,I>::internal_preprocess_submit()
{
    data->d_perm_to_sort_back_sc.alloc();

    gpu::mgm::copy_submit(q, data->h_perm_to_sort_back_sc, data->d_perm_to_sort_back_sc);

    data->op_X_sp2dn->set_ws_persistent(ator_ws_persistent->alloc(data->op_X_sp2dn->get_wss_persistent()));
    data->op_X_sp2dn->preprocess_submit(ator_ws_tmp_overlap->alloc(data->op_X_sp2dn->get_wss_tmp_preprocess()));

    data->op_d_trsm.set_ws_persistent(ator_ws_persistent->alloc(data->op_d_trsm.get_wss_persistent()));
    data->op_d_trsm.preprocess_submit(ator_ws_tmp_overlap->alloc(data->op_d_trsm.get_wss_tmp_preprocess()));
}



template<typename T, typename I>
void sc_hcsx_ddny_tria<T,I>::internal_perform_1_submit()
{
    stacktimer::push("numerical_factorization");
    data->h_A11_solver.numericalFactorization();
    stacktimer::pop();
}



template<typename T, typename I>
void sc_hcsx_ddny_tria<T,I>::internal_perform_2_submit()
{
    stacktimer::push("extract_factors");
    if(data->solver_factor_uplo == 'U') {
        Matrix_CSR<T,I,gpu::mgm::Ah> h_factor_old = MatrixCsxView_new<T,I>::template to_old<gpu::mgm::Ah>(data->h_factor_U_row);
        data->h_A11_solver.getFactorU(h_factor_old, false, true);
    }
    if(data->solver_factor_uplo == 'L') {
        Matrix_CSR<T,I,gpu::mgm::Ah> h_factor_old = MatrixCsxView_new<T,I>::template to_old<gpu::mgm::Ah>(data->h_factor_L_row);
        data->h_A11_solver.getFactorL(h_factor_old, false, true);
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
void sc_hcsx_ddny_tria<T,I>::internal_solve_A11(VectorDenseView_new<T> & rhs, VectorDenseView_new<T> & sol)
{
    Vector_Dense<T,I> rhs_old = VectorDenseView_new<T>::template to_old<I,cpu_allocator>(rhs);
    Vector_Dense<T,I> sol_old = VectorDenseView_new<T>::template to_old<I,cpu_allocator>(sol);

    data->h_A11_solver.solve(rhs_old, sol_old);
}



#define INSTANTIATE_T_I(T,I) \
template class sc_hcsx_ddny_tria<T,I>;

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
