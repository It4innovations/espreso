
#include "math/operations/sc_csx_dny.tria.h"

#include "esinfo/meshinfo.h"
#include "basis/containers/allocators.h"
#include "math/primitives/matrix_csr.h"

#include "math/operations/trsm_csx_dny_tri.h"
#include "math/operations/herk_dnx_dny_tri.h"
#include "math/operations/convert_csx_csy_map.h"
#include "math/operations/pivots_trails_csx.h"
#include "math/operations/sorting_permutation.h"
#include "math/operations/convert_csx_dny.h"
#include "math/operations/complete_dnx.h"
#include "math/operations/permute_dnx_dnx.h"
#include "math/operations/permute_csx_csx.h"
#include "math/operations/copy_dnx.h"
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
    };
    config cfg;
    Matrix_CSR<T,I> A11_old;
    DirectSparseSolver<T,I> A11_solver;
    bool need_reorder_factor_L2U = false;
    bool need_reorder_factor_U2L = false;
    char solver_factor_uplo = '_';
    MatrixCsxData_new<T,I> factor_L_row;
    MatrixCsxData_new<T,I> factor_U_row;
    MatrixCsxView_new<T,I> L_row;
    MatrixCsxView_new<T,I> L_col;
    MatrixCsxView_new<T,I> U_row;
    MatrixCsxView_new<T,I> U_col;
    MatrixCsxView_new<T,I> * L_to_use = nullptr;
    MatrixCsxData_new<T,I> X_sp;
    MatrixDenseData_new<T> X_dn;
    MatrixDenseData_new<T> sc_tmp1;
    MatrixDenseData_new<T> sc_tmp2;
    PermutationData_new<I> perm_to_sort_A12_cols;
    PermutationView_new<I> perm_to_sort_back_sc;
    PermutationData_new<I> perm_fillreduce;
    convert_csx_csy_map<T,I> op_L2U;
    convert_csx_csy_map<T,I> op_U2L;
    trsm_csx_dny_tri<T,I> op_trsm;
    herk_dnx_dny_tri<T,I> op_herk;
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
    stacktimer::push("sc_csx_dny_tria::preprocess");

    if(!is_matrix_hermitian) eslog::error("dont support non-hermitian systems yet\n");
    if(A12 == nullptr) eslog::error("A12 has to be set for now\n");

    data = std::make_unique<sc_csx_dny_tria_data<T,I>>();

    populate_config<T,I>(data->cfg, size_A11);

    if(!DirectSparseSolver<T,I>::provideFactors()) eslog::error("wrong sparse solver, must provide factors\n");
    data->solver_factor_uplo = '_';
    if(DirectSparseSolver<T,I>::factorsSymmetry() == Solver_Factors::HERMITIAN_LOWER) data->solver_factor_uplo = 'L';
    if(DirectSparseSolver<T,I>::factorsSymmetry() == Solver_Factors::HERMITIAN_UPPER) data->solver_factor_uplo = 'U';
    if(data->solver_factor_uplo == '_') eslog::error("wrong sparse solver, must be symmetric\n");

    if(called_set_matrix == '1') {
        eslog::error("currrently dont support single large matrix on input\n");
        // if hermitian then be carefull about A12 A21 which of them is full
    }

    if(A11->order != 'R') {
        eslog::error("only support csr matrices now\n");
    }

    if(A22 != nullptr) {
        eslog::error("for now I assume A22=O\n");
    }

    stacktimer::push("symbolic_factorization");
    data->A11_old = MatrixCsxView_new<T,I>::template to_old<cpu_allocator>(*A11);
    data->A11_solver.commit(data->A11_old);
    data->A11_solver.symbolicFactorization();
    stacktimer::pop();

    I factor_nnz = data->A11_solver.getFactorNnz();

    data->factor_U_row.set(size_A11, size_A11, factor_nnz, 'R', AllocatorCPU_new::get_singleton());
    data->factor_U_row.prop.uplo = 'U';
    data->factor_U_row.prop.diag = 'N';

    data->factor_L_row.set(size_A11, size_A11, factor_nnz, 'R', AllocatorCPU_new::get_singleton());
    data->factor_L_row.prop.uplo = 'L';
    data->factor_L_row.prop.diag = 'N';

    data->need_reorder_factor_L2U = ((data->solver_factor_uplo == 'L') && (data->cfg.order_L == 'C'));
    data->need_reorder_factor_U2L = ((data->solver_factor_uplo == 'U') && (data->cfg.order_L == 'R'));

    if(data->solver_factor_uplo == 'U' || data->need_reorder_factor_L2U) {
        data->factor_U_row.alloc();
    }
    if(data->solver_factor_uplo == 'L' || data->need_reorder_factor_U2L) {
        data->factor_L_row.alloc();
    }

    stacktimer::push("extract_factors");
    if(data->solver_factor_uplo == 'U') {
        Matrix_CSR<T,I> factor_old_U = MatrixCsxView_new<T,I>::template to_old<cpu_allocator>(data->factor_U_row);
        data->A11_solver.getFactorU(factor_old_U, true, false);
    }
    if(data->solver_factor_uplo == 'L') {
        Matrix_CSR<T,I> factor_old_L = MatrixCsxView_new<T,I>::template to_old<cpu_allocator>(data->factor_L_row);
        data->A11_solver.getFactorL(factor_old_L, true, false);
    }
    stacktimer::pop();

    data->L_row = data->factor_L_row;
    data->L_col = data->factor_U_row.get_transposed_reordered_view();
    data->U_row = data->factor_U_row;
    data->U_col = data->factor_L_row.get_transposed_reordered_view();

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
            Permutation<I> perm_fillreduce_old = PermutationView_new<I>::template to_old<cpu_allocator>(data->perm_fillreduce);
            data->A11_solver.getPermutation(perm_fillreduce_old);

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

    data->sc_tmp1.set(sc->nrows, sc->ncols, sc->order, AllocatorCPU_new::get_singleton());
    data->sc_tmp1.prop.uplo = sc->prop.uplo;

    data->sc_tmp2.set(sc->nrows, sc->ncols, sc->order, AllocatorCPU_new::get_singleton());
    data->sc_tmp2.prop.uplo = sc->prop.uplo;

    data->op_trsm.set_config(data->cfg.cfg_trsm);
    data->op_trsm.set_L(data->L_to_use);
    data->op_trsm.set_X(&data->X_dn);
    data->op_trsm.calc_X_pattern(data->X_sp);
    data->op_trsm.preprocess();

    data->op_herk.set_config(data->cfg.cfg_herk);
    data->op_herk.set_matrix_A(&data->X_dn);
    data->op_herk.set_matrix_C(&data->sc_tmp1);
    data->op_herk.set_coefficients(-alpha, 0 * alpha);  // beta=0, because I assume A22=0
    data->op_herk.set_mode(blas::herk_mode::AhA);
    data->op_herk.calc_A_pattern(data->X_sp);
    data->op_herk.preprocess();

    stacktimer::pop();
}



template<typename T, typename I>
void sc_csx_dny_tria<T,I>::internal_perform_1()
{
    stacktimer::push("sc_csx_dny_tria::perform_1");

    stacktimer::push("numerical_factorization");
    data->A11_solver.numericalFactorization();
    stacktimer::pop();

    stacktimer::pop();
}



template<typename T, typename I>
void sc_csx_dny_tria<T,I>::internal_perform_2()
{
    stacktimer::push("sc_csx_dny_tria::perform_2");

    stacktimer::push("extract_factors");
    if(data->solver_factor_uplo == 'U') {
        Matrix_CSR<T,I,gpu::mgm::Ah> factor_old = MatrixCsxView_new<T,I>::template to_old<gpu::mgm::Ah>(data->factor_U_row);
        data->A11_solver.getFactorU(factor_old, false, true);
    }
    if(data->solver_factor_uplo == 'L') {
        Matrix_CSR<T,I,gpu::mgm::Ah> factor_old = MatrixCsxView_new<T,I>::template to_old<gpu::mgm::Ah>(data->factor_L_row);
        data->A11_solver.getFactorL(factor_old, false, true);
    }
    stacktimer::pop();

    if(data->need_reorder_factor_L2U) {
        data->op_L2U.perform_values();
    }
    if(data->need_reorder_factor_U2L) {
        data->op_U2L.perform_values();
    }

    data->X_dn.alloc();
    data->sc_tmp1.alloc();
    data->sc_tmp2.alloc();

    convert_csx_dny<T,I>::do_all(&data->X_sp, &data->X_dn);

    data->op_trsm.perform();

    data->op_herk.perform();

    complete_dnx<T>::do_all(&data->sc_tmp1, sc->prop.uplo, true);

    permute_dnx_dnx<T,I>::do_all(&data->sc_tmp1, &data->sc_tmp2, &data->perm_to_sort_back_sc, &data->perm_to_sort_back_sc);

    copy_dnx<T>::do_all(&data->sc_tmp2, sc, false);

    data->X_dn.free();
    data->sc_tmp1.free();
    data->sc_tmp2.free();

    stacktimer::pop();
}



template<typename T, typename I>
void sc_csx_dny_tria<T,I>::internal_solve_A11(VectorDenseView_new<T> & rhs, VectorDenseView_new<T> & sol)
{
    Vector_Dense<T,I> rhs_old = VectorDenseView_new<T>::template to_old<I,cpu_allocator>(rhs);
    Vector_Dense<T,I> sol_old = VectorDenseView_new<T>::template to_old<I,cpu_allocator>(sol);

    data->A11_solver.solve(rhs_old, sol_old);
}



#define INSTANTIATE_T_I(T,I) \
template class sc_csx_dny_tria<T,I>;

    #define INSTANTIATE_T(T) \
    INSTANTIATE_T_I(T, int32_t) \
    /* INSTANTIATE_T_I(T, int64_t) */

        /* INSTANTIATE_T(float) */ \
        INSTANTIATE_T(double) \
        /* INSTANTIATE_T(std::complex<float>) */ \
        INSTANTIATE_T(std::complex<double>)

    #undef INSTANTIATE_T
#undef INSTANTIATE_T_I



}
}
}
