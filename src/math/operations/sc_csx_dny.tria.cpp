
#include "math/operations/sc_csx_dny.tria.h"

#include "esinfo/meshinfo.h"
#include "basis/containers/allocators.h"
#include "math/primitives/matrix_csr.h"

#include "math/primitives_new/allocator_new.h"
#include "math/operations/trsm_csx_dny_tri.h"
#include "math/operations/herk_dnx_dny_tri.h"
#include "math/operations/convert_csx_csy_map.h"
#include "math/operations/pivots_trails_csx.h"
#include "math/operations/sorting_permutation.h"
#include "math/operations/convert_csx_dny.h"
#include "math/operations/complete_dnx_dnx.h"
#include "math/operations/permute_dnx_dnx.h"
#include "math/operations/permute_csx_csx.h"
#include "math/operations/copy_dnx.h"
#include "math/operations/quadrisect_csx_csy.h"
#include "math/operations/lincomb_dnx_csy.h"
#include "math/operations/solver_csx.h"
#include "basis/utilities/stacktimer.h"



namespace espreso {
namespace math {
namespace operations {



template<typename T, typename I>
struct sc_csx_dny_tria_data
{
    struct config
    {
        typename trsm_csx_dny_tri<T,I>::config cfg_trsm;
        typename herk_dnx_dny_tri<T,I>::config cfg_herk;
        char order_X = '_';
        char order_L = '_';
        typename solver_csx<T,I>::implementation_selector op_a11solver_is = solver_csx<T,I>::implementation_selector::autoselect;
    };
    config cfg;
    bool need_reorder_factor_L2U = false;
    bool need_reorder_factor_U2L = false;
    char solver_factor_uplo = '_';
    MatrixCsxData_new<T,I> L_data;
    MatrixCsxData_new<T,I> U_data;
    MatrixCsxView_new<T,I> L_row;
    MatrixCsxView_new<T,I> L_col;
    MatrixCsxView_new<T,I> U_row;
    MatrixCsxView_new<T,I> U_col;
    MatrixCsxView_new<T,I> * L_to_use = nullptr;
    MatrixCsxData_new<T,I> X_sp;
    MatrixDenseData_new<T> X_dn;
    PermutationData_new<I> perm_to_sort_A12_cols;
    PermutationView_new<I> perm_to_sort_back_sc;
    PermutationData_new<I> perm_fillreduce;
    std::unique_ptr<solver_csx<T,I>> op_A11_solver;
    convert_csx_csy_map<T,I> op_L2U;
    convert_csx_csy_map<T,I> op_U2L;
    trsm_csx_dny_tri<T,I> op_trsm;
    herk_dnx_dny_tri<T,I> op_herk;
    permute_dnx_dnx<T,I> op_permute_sc;
    lincomb_dnx_csy<T,I> op_lincomb_final;
    MatrixCsxView_new<T,I> A22_rt;
    quadrisect_csx_csy<T,I> op_split;
    MatrixCsxData_new<T,I> sub_A11;
    MatrixCsxData_new<T,I> sub_A12;
    MatrixCsxData_new<T,I> sub_A21;
    MatrixCsxData_new<T,I> sub_A22;
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
static void populate_config_from_env(typename sc_csx_dny_tria_data<T,I>::config & cfg)
{
    cfg.cfg_trsm.partition.parameter = 0;
    cfg.cfg_herk.partition_parameter = 0;

    set_by_env(cfg.order_X,                                   "ESPRESO_CONFIG_OP_SC_CSX_DNY_TRIA_order_X");
    set_by_env(cfg.order_L,                                   "ESPRESO_CONFIG_OP_SC_CSX_DNY_TRIA_order_L");
    set_by_env(cfg.cfg_trsm.strategy,                         "ESPRESO_CONFIG_OP_SC_CSX_DNY_TRIA_trsm_strategy");
    set_by_env(cfg.cfg_trsm.partition.algorithm,              "ESPRESO_CONFIG_OP_SC_CSX_DNY_TRIA_trsm_partition_algorithm");
    set_by_env(cfg.cfg_trsm.partition.parameter,              "ESPRESO_CONFIG_OP_SC_CSX_DNY_TRIA_trsm_partition_parameter");
    set_by_env(cfg.cfg_trsm.splitrhs.factor_order_sp,         "ESPRESO_CONFIG_OP_SC_CSX_DNY_TRIA_trsm_splitrhs_factor_order_sp");
    set_by_env(cfg.cfg_trsm.splitrhs.factor_order_dn,         "ESPRESO_CONFIG_OP_SC_CSX_DNY_TRIA_trsm_splitrhs_factor_order_dn");
    set_by_env(cfg.cfg_trsm.splitrhs.spdn_criteria,           "ESPRESO_CONFIG_OP_SC_CSX_DNY_TRIA_trsm_splitrhs_spdn_criteria");
    set_by_env(cfg.cfg_trsm.splitrhs.spdn_param,              "ESPRESO_CONFIG_OP_SC_CSX_DNY_TRIA_trsm_splitrhs_spdn_param");
    set_by_env(cfg.cfg_trsm.splitfactor.trsm_factor_spdn,     "ESPRESO_CONFIG_OP_SC_CSX_DNY_TRIA_trsm_splitfactor_trsm_factor_spdn");
    set_by_env(cfg.cfg_trsm.splitfactor.trsm_factor_order,    "ESPRESO_CONFIG_OP_SC_CSX_DNY_TRIA_trsm_splitfactor_trsm_factor_order");
    set_by_env(cfg.cfg_trsm.splitfactor.gemm_factor_order_sp, "ESPRESO_CONFIG_OP_SC_CSX_DNY_TRIA_trsm_splitfactor_gemm_factor_order_sp");
    set_by_env(cfg.cfg_trsm.splitfactor.gemm_factor_order_dn, "ESPRESO_CONFIG_OP_SC_CSX_DNY_TRIA_trsm_splitfactor_gemm_factor_order_dn");
    set_by_env(cfg.cfg_trsm.splitfactor.gemm_factor_prune,    "ESPRESO_CONFIG_OP_SC_CSX_DNY_TRIA_trsm_splitfactor_gemm_factor_prune");
    set_by_env(cfg.cfg_trsm.splitfactor.gemm_spdn_criteria,   "ESPRESO_CONFIG_OP_SC_CSX_DNY_TRIA_trsm_splitfactor_gemm_spdn_criteria");
    set_by_env(cfg.cfg_trsm.splitfactor.gemm_spdn_param,      "ESPRESO_CONFIG_OP_SC_CSX_DNY_TRIA_trsm_splitfactor_gemm_spdn_param");
    set_by_env(cfg.cfg_herk.strategy,                         "ESPRESO_CONFIG_OP_SC_CSX_DNY_TRIA_herk_strategy");
    set_by_env(cfg.cfg_herk.partition_algorithm,              "ESPRESO_CONFIG_OP_SC_CSX_DNY_TRIA_herk_partition_algorithm");
    set_by_env(cfg.cfg_herk.partition_parameter,              "ESPRESO_CONFIG_OP_SC_CSX_DNY_TRIA_herk_partition_parameter");
}

template<typename T, typename I>
static void populate_config_replace_defaults(typename sc_csx_dny_tria_data<T,I>::config & cfg, size_t ndofs_domain)
{
    replace_if_default(cfg.cfg_trsm.strategy, 'F');
    bool is_in_between = ((ndofs_domain > 1000) && (ndofs_domain < 16000));
    replace_if_default(cfg.cfg_herk.strategy, is_in_between ? 'Q' : 'T');

    if(info::mesh->dimension == 2 && cfg.cfg_herk.strategy == 'Q') replace_if_zero(cfg.cfg_herk.partition_parameter, -200);
    if(info::mesh->dimension == 2 && cfg.cfg_herk.strategy == 'T') replace_if_zero(cfg.cfg_herk.partition_parameter, -200);
    if(info::mesh->dimension == 3 && cfg.cfg_herk.strategy == 'Q') replace_if_zero(cfg.cfg_herk.partition_parameter, 50);
    if(info::mesh->dimension == 3 && cfg.cfg_herk.strategy == 'T') replace_if_zero(cfg.cfg_herk.partition_parameter, 10);
    replace_if_default(cfg.cfg_herk.partition_algorithm, 'U');

    if(info::mesh->dimension == 2 && cfg.cfg_trsm.strategy == 'F') replace_if_zero(cfg.cfg_trsm.partition.parameter, -200);
    if(info::mesh->dimension == 2 && cfg.cfg_trsm.strategy == 'R') replace_if_zero(cfg.cfg_trsm.partition.parameter, -100);
    if(info::mesh->dimension == 3 && cfg.cfg_trsm.strategy == 'F') replace_if_zero(cfg.cfg_trsm.partition.parameter, -200);
    if(info::mesh->dimension == 3 && cfg.cfg_trsm.strategy == 'R') replace_if_zero(cfg.cfg_trsm.partition.parameter, -100);
    replace_if_default(cfg.cfg_trsm.partition.algorithm, 'U');

    replace_if_default(cfg.cfg_trsm.splitrhs.spdn_criteria, 'S');
    replace_if_default(cfg.cfg_trsm.splitfactor.trsm_factor_spdn, (info::mesh->dimension == 2) ? 'S' : 'D');
    replace_if_default(cfg.cfg_trsm.splitfactor.gemm_factor_prune, 'R');
    replace_if_default(cfg.cfg_trsm.splitfactor.gemm_spdn_criteria, (info::mesh->dimension == 3 && cfg.cfg_trsm.splitfactor.gemm_factor_prune == 'R') ? 'D' : 'S');

    replace_if_default(cfg.order_X, (info::mesh->dimension == 2) ? 'C' : 'R');
    replace_if_default(cfg.cfg_trsm.splitrhs.factor_order_sp, 'R');
    replace_if_default(cfg.cfg_trsm.splitrhs.factor_order_dn, 'C');
    replace_if_default(cfg.cfg_trsm.splitfactor.trsm_factor_order, 'R');
    replace_if_default(cfg.cfg_trsm.splitfactor.gemm_factor_order_sp, 'R');
    replace_if_default(cfg.cfg_trsm.splitfactor.gemm_factor_order_dn, 'R');

    replace_if_default(cfg.order_L, 'C');
}

template<typename T, typename I>
static void populate_config(typename sc_csx_dny_tria_data<T,I>::config & cfg, size_t ndofs_domain)
{
    populate_config_from_env<T,I>(cfg);

    populate_config_replace_defaults<T,I>(cfg, ndofs_domain);
}



template<typename T, typename I>
sc_csx_dny_tria<T,I>::sc_csx_dny_tria() = default;



template<typename T, typename I>
sc_csx_dny_tria<T,I>::~sc_csx_dny_tria() = default;



template<typename T, typename I>
void sc_csx_dny_tria<T,I>::internal_preprocess()
{
    if(!is_matrix_hermitian) eslog::error("dont support non-hermitian systems yet\n");

    data = std::make_unique<sc_csx_dny_tria_data<T,I>>();

    populate_config<T,I>(data->cfg, size_A11);

    if(called_set_matrix == '1') {
        data->op_split.set_matrix_src(A);
        data->op_split.set_matrices_dst(&data->sub_A11, &data->sub_A12, &data->sub_A21, &data->sub_A22);
        data->op_split.set_bounds(size_A11, size_A11);
        data->op_split.setup();

        data->sub_A11.set(size_A11, size_A11, data->op_split.get_output_matrix_11_nnz(), A->order, AllocatorCPU_new::get_singleton());
        data->sub_A12.set(size_A11, size_sc,  data->op_split.get_output_matrix_12_nnz(), A->order, AllocatorCPU_new::get_singleton());
        data->sub_A21.set(size_sc,  size_A11, data->op_split.get_output_matrix_21_nnz(), A->order, AllocatorCPU_new::get_singleton());
        data->sub_A22.set(size_sc,  size_sc,  data->op_split.get_output_matrix_22_nnz(), A->order, AllocatorCPU_new::get_singleton());

        data->sub_A11.prop = A->prop;
        data->sub_A22.prop = A->prop;

        data->sub_A11.alloc();
        data->sub_A12.alloc();
        data->sub_A21.alloc();
        data->sub_A22.alloc();

        data->op_split.perform();

        A11 = &data->sub_A11;
        if(A->prop.uplo != 'L') A12 = &data->sub_A12;
        if(A->prop.uplo != 'U') A21 = &data->sub_A21;
        A22 = &data->sub_A22;
    }

    data->solver_factor_uplo = A11->prop.uplo;

    if(A12 == nullptr) eslog::error("A12 has to be set for now\n");

    if(A11->order != 'R') {
        eslog::error("only support csr matrices now\n");
    }

    data->op_A11_solver = solver_csx<T,I>::make(data->cfg.op_a11solver_is, A11, true, need_solve_A11);
    data->op_A11_solver->set_matrix_A(A11);
    data->op_A11_solver->set_needs(true, need_solve_A11);
    data->op_A11_solver->factorize_symbolic();

    size_t factor_nnz = data->op_A11_solver->get_factor_nnz();

    data->U_data.set(size_A11, size_A11, factor_nnz, A11->order, AllocatorCPU_new::get_singleton());
    data->U_data.prop.uplo = 'U';
    data->U_data.prop.diag = 'N';

    data->L_data.set(size_A11, size_A11, factor_nnz, A11->order, AllocatorCPU_new::get_singleton());
    data->L_data.prop.uplo = 'L';
    data->L_data.prop.diag = 'N';

    data->need_reorder_factor_L2U = ((data->solver_factor_uplo == 'L') && (data->cfg.order_L == 'C'));
    data->need_reorder_factor_U2L = ((data->solver_factor_uplo == 'U') && (data->cfg.order_L == 'R'));

    if(data->solver_factor_uplo == 'U' || data->need_reorder_factor_L2U) {
        data->U_data.alloc();
    }
    if(data->solver_factor_uplo == 'L' || data->need_reorder_factor_U2L) {
        data->L_data.alloc();
    }

    if(A11->prop.uplo == 'L') {
        data->op_A11_solver->get_factor_L(data->L_data, true, false);
    }
    if(A11->prop.uplo == 'U') {
        data->op_A11_solver->get_factor_U(data->U_data, true, false);
    }

    if(A11->order == 'R') {
        data->L_row = data->L_data;
        data->U_row = data->U_data;
        data->L_col = data->U_row.get_transposed_reordered_view();
        data->U_col = data->L_row.get_transposed_reordered_view();
    }
    if(A11->order == 'C') {
        data->L_col = data->L_data;
        data->U_col = data->U_data;
        data->L_row = data->U_col.get_transposed_reordered_view();
        data->U_row = data->L_col.get_transposed_reordered_view();
    }

    if(data->need_reorder_factor_L2U) {
        data->op_L2U.set_matrix_src(&data->L_row);
        data->op_L2U.set_matrix_dst(&data->L_col);
        data->op_L2U.perform_pattern();
    }
    if(data->need_reorder_factor_U2L) {
        data->op_U2L.set_matrix_src(&data->U_row);
        data->op_U2L.set_matrix_dst(&data->U_col);
        data->op_U2L.perform_pattern();
    }

    if(data->cfg.order_L == 'R') data->L_to_use = &data->L_row;
    if(data->cfg.order_L == 'C') data->L_to_use = &data->L_col;

    data->X_sp.set(A12->nrows, A12->ncols, A12->nnz, A12->order, AllocatorCPU_new::get_singleton());
    data->X_sp.alloc();

    {
        MatrixCsxData_new<T,I> X_sp_tmp;
        X_sp_tmp.set(data->X_sp.nrows, data->X_sp.ncols, data->X_sp.nnz, data->X_sp.order, AllocatorCPU_new::get_singleton());
        X_sp_tmp.alloc();

        {
            data->perm_fillreduce.set(A12->nrows, AllocatorCPU_new::get_singleton());
            data->perm_fillreduce.alloc();
            data->op_A11_solver->get_permutation(data->perm_fillreduce);

            math::operations::permute_csx_csx<T,I>::do_all(A12, &X_sp_tmp, &data->perm_fillreduce, nullptr);
        }

        {
            VectorDenseData_new<I> colpivots;
            colpivots.set(A12->ncols, AllocatorCPU_new::get_singleton());
            colpivots.alloc();

            math::operations::pivots_trails_csx<T,I>::do_all(&X_sp_tmp, &colpivots, 'C', 'P', '_');

            data->perm_to_sort_A12_cols.set(A12->ncols, AllocatorCPU_new::get_singleton());
            data->perm_to_sort_A12_cols.alloc();

            data->perm_to_sort_back_sc = data->perm_to_sort_A12_cols.get_inverse_view();

            math::operations::sorting_permutation<I,I>::do_all(&colpivots, &data->perm_to_sort_A12_cols);
        }

        math::operations::permute_csx_csx<T,I>::do_all(&X_sp_tmp, &data->X_sp, nullptr, &data->perm_to_sort_A12_cols);
    }

    data->X_dn.set(data->X_sp.nrows, data->X_sp.ncols, data->cfg.order_X, AllocatorCPU_new::get_singleton());

    data->op_trsm.set_config(data->cfg.cfg_trsm);
    data->op_trsm.set_L(data->L_to_use);
    data->op_trsm.set_X(&data->X_dn);
    data->op_trsm.calc_X_pattern(data->X_sp);
    data->op_trsm.preprocess();

    data->op_herk.set_config(data->cfg.cfg_herk);
    data->op_herk.set_matrix_A(&data->X_dn);
    data->op_herk.set_matrix_C(sc);
    data->op_herk.set_coefficients(-alpha, 0);
    data->op_herk.set_mode(blas::herk_mode::AhA);
    data->op_herk.calc_A_pattern(data->X_sp);
    data->op_herk.preprocess();

    data->op_permute_sc.set_matrix_src(sc);
    data->op_permute_sc.set_matrix_dst(sc);
    data->op_permute_sc.set_perm_vector_rows(&data->perm_to_sort_back_sc);
    data->op_permute_sc.set_perm_vector_cols(&data->perm_to_sort_back_sc);

    if(A22 != nullptr) {
        data->A22_rt = A22->get_transposed_reordered_view();

        data->op_lincomb_final.set_matrix_X(sc);
        data->op_lincomb_final.set_matrix_A(sc);
        if(A22->prop.uplo == sc->prop.uplo) data->op_lincomb_final.set_matrix_B(A22);
        else data->op_lincomb_final.set_matrix_B(&data->A22_rt);
        data->op_lincomb_final.set_coefficients(T{1}, alpha);
    }
}



template<typename T, typename I>
void sc_csx_dny_tria<T,I>::internal_perform_1()
{
    if(called_set_matrix == '1') {
        data->op_split.perform();
    }

    data->op_A11_solver->factorize_numeric();
}



template<typename T, typename I>
void sc_csx_dny_tria<T,I>::internal_perform_2()
{
    if(data->solver_factor_uplo == 'U') {
        data->op_A11_solver->get_factor_U(data->U_data, false, true);
    }
    if(data->solver_factor_uplo == 'L') {
        data->op_A11_solver->get_factor_L(data->L_data, false, true);
    }

    if(data->need_reorder_factor_L2U) {
        data->op_L2U.perform_values();
    }
    if(data->need_reorder_factor_U2L) {
        data->op_U2L.perform_values();
    }

    data->X_dn.alloc();

    permute_csx_csx<T,I>::do_all(A12, &data->X_sp, &data->perm_fillreduce, &data->perm_to_sort_A12_cols);

    convert_csx_dny<T,I>::do_all(&data->X_sp, &data->X_dn);

    data->op_trsm.perform();

    data->op_herk.perform();

    data->op_permute_sc.perform();

    data->X_dn.free();

    if(A22 != nullptr) {
        data->A22_rt = A22->get_transposed_reordered_view();

        data->op_lincomb_final.perform();
    }
}



template<typename T, typename I>
void sc_csx_dny_tria<T,I>::internal_solve_A11(VectorDenseView_new<T> & rhs, VectorDenseView_new<T> & sol)
{
    data->op_A11_solver->solve(rhs, sol);
}



#define INSTANTIATE_T_I(T,I) \
template class sc_csx_dny_tria<T,I>;

    #define INSTANTIATE_T(T) \
    INSTANTIATE_T_I(T, int32_t) \
    /* INSTANTIATE_T_I(T, int64_t) */

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
