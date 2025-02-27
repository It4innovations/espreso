#include "math/wrappers/math.sc_solver.h"

#ifdef ESPRESO_USE_WRAPPER_SCSOLVER_SCSOLVERTRIANGULAR

#include <cstdlib>

#include "math/primitives_new/matrix_csx_data_new.h"
#include "math/primitives_new/matrix_dense_data_new.h"
#include "math/primitives_new/permutation_data_new.h"
#include "math/operations/convert_csx_csy_map.h"
#include "math/operations/trsm_csx_dny_tri.h"
#include "math/operations/herk_dnx_dny_tri.h"
#include "math/operations/permute_csx_csx.h"
#include "math/operations/permute_dnx_dnx.h"
#include "math/operations/pivots_trails_csx.h"
#include "math/operations/sorting_permutation.h"
#include "math/operations/convert_csx_dny.h"
#include "math/operations/complete_dnx.h"
#include "math/operations/lincomb_matrix_dnx.h"



namespace espreso {



template<typename T, typename I>
static void extract_submatrices(const Matrix_CSR<T,I> & A, Matrix_CSR<T,I> & A11, Matrix_CSR<T,I> & A12, Matrix_CSR<T,I> & A21, Matrix_CSR<T,I> & A22, I sc_size)
{
    I other_size = A.nrows - sc_size;

    math::spblas::submatrix(A, A11, 0,          other_size, 0,          other_size);
    math::spblas::submatrix(A, A12, 0,          other_size, other_size, A.nrows);
    math::spblas::submatrix(A, A21, other_size, A.nrows,    0,          other_size);
    math::spblas::submatrix(A, A22, other_size, A.nrows,    other_size, A.nrows);
}



template<typename T, typename I>
struct Schur_Complement_Solver_External_Representation
{
    const Matrix_CSR<T,I> * A_whole = nullptr;
    const Matrix_CSR<T,I> * A11 = nullptr;
    const Matrix_CSR<T,I> * A12 = nullptr;
    const Matrix_CSR<T,I> * A21 = nullptr;
    const Matrix_CSR<T,I> * A22 = nullptr;
    Matrix_CSR<T,I> A11_my;
    Matrix_CSR<T,I> A12_my;
    Matrix_CSR<T,I> A21_my;
    Matrix_CSR<T,I> A22_my;
    MatrixCsxView_new<T,I> A22_new;
    DirectSparseSolver<T,I> A11_solver;
    Solver_Factors A11_solver_factors;
    bool A11_solver_factors_get_L;
    bool A11_solver_factors_get_U;
    bool system_is_hermitian;
    Matrix_CSR<T,I> L11;
    Matrix_CSR<T,I> U11;
    Permutation<I> perm;
    PermutationView_new<I> perm_new;
    PermutationData_new<I> perm_to_sort;
    MatrixCsxView_new<T,I> L_new;
    MatrixCsxView_new<T,I> U_new;
    MatrixCsxData_new<T,I> L_new_reordered;
    MatrixCsxView_new<T,I> * L_new_to_use;
    MatrixCsxView_new<T,I> B_right_input;
    MatrixCsxData_new<T,I> B_right; // already permuted and sorted
    MatrixDenseData_new<T> X;
    math::operations::convert_csx_csy_map<T,I> op_reorder_L;
    math::operations::trsm_csx_dny_tri<T,I> op_trsm_tri_L;
    math::operations::herk_dnx_dny_tri<T,I> op_herk_tri;
    I sc_size = 0;
    I other_size = 0;
    int stage = 0;
    struct config {
        typename math::operations::trsm_csx_dny_tri<T,I>::config trsm_tri_cfg;
        typename math::operations::herk_dnx_dny_tri<T,I>::config herk_tri_cfg;
        char order_X = '_';
        char order_L = '_';
    } cfg;
};



template<typename T>
static void set_by_env(T var, const char * env_var)
{
    const char * env_val = std::getenv(env_var);
    if(env_val != nullptr) {
        if constexpr(std::is_same_v<T,char>) var = *env_val;
        if constexpr(std::is_same_v<T,int>) var = atoi(env_val);
        if constexpr(std::is_same_v<T,double>) var = atof(env_val);
    }
}



template<typename T, typename I>
static void init_config(typename Schur_Complement_Solver_External_Representation<T,I>::config & cfg)
{
    set_by_env(cfg.order_X, "ESPRESO_SCSOLVERTRIANGULAR_CONFIG_order_X");
    set_by_env(cfg.order_L, "ESPRESO_SCSOLVERTRIANGULAR_CONFIG_order_L");

    set_by_env(cfg.trsm_tri_cfg.strategy, "ESPRESO_SCSOLVERTRIANGULAR_CONFIG_trsm_strategy");
    set_by_env(cfg.trsm_tri_cfg.partition.algorithm, "ESPRESO_SCSOLVERTRIANGULAR_CONFIG_trsm_partition_algorithm");
    set_by_env(cfg.trsm_tri_cfg.partition.parameter, "ESPRESO_SCSOLVERTRIANGULAR_CONFIG_trsm_partition_parameter");
    set_by_env(cfg.trsm_tri_cfg.splitrhs.factor_order_sp, "ESPRESO_SCSOLVERTRIANGULAR_CONFIG_trsm_splitrhs_factor_order_sp");
    set_by_env(cfg.trsm_tri_cfg.splitrhs.factor_order_dn, "ESPRESO_SCSOLVERTRIANGULAR_CONFIG_trsm_splitrhs_factor_order_dn");
    set_by_env(cfg.trsm_tri_cfg.splitrhs.spdn_criteria, "ESPRESO_SCSOLVERTRIANGULAR_CONFIG_trsm_splitrhs_spdn_criteria");
    set_by_env(cfg.trsm_tri_cfg.splitrhs.spdn_param, "ESPRESO_SCSOLVERTRIANGULAR_CONFIG_trsm_splitrhs_spdn_param");
    set_by_env(cfg.trsm_tri_cfg.splitfactor.trsm_factor_spdn, "ESPRESO_SCSOLVERTRIANGULAR_CONFIG_trsm_splitfactor_trsm_factor_spdn");
    set_by_env(cfg.trsm_tri_cfg.splitfactor.trsm_factor_order, "ESPRESO_SCSOLVERTRIANGULAR_CONFIG_trsm_splitfactor_trsm_factor_order");
    set_by_env(cfg.trsm_tri_cfg.splitfactor.gemm_factor_prune, "ESPRESO_SCSOLVERTRIANGULAR_CONFIG_trsm_splitfactor_gemm_factor_prune");
    set_by_env(cfg.trsm_tri_cfg.splitfactor.gemm_factor_order_sp, "ESPRESO_SCSOLVERTRIANGULAR_CONFIG_trsm_splitfactor_gemm_factor_order_sp");
    set_by_env(cfg.trsm_tri_cfg.splitfactor.gemm_factor_order_dn, "ESPRESO_SCSOLVERTRIANGULAR_CONFIG_trsm_splitfactor_gemm_factor_order_dn");
    set_by_env(cfg.trsm_tri_cfg.splitfactor.gemm_spdn_criteria, "ESPRESO_SCSOLVERTRIANGULAR_CONFIG_trsm_splitfactor_gemm_spdn_criteria");
    set_by_env(cfg.trsm_tri_cfg.splitfactor.gemm_spdn_param, "ESPRESO_SCSOLVERTRIANGULAR_CONFIG_trsm_splitfactor_gemm_spdn_param");

    set_by_env(cfg.herk_tri_cfg.strategy, "ESPRESO_SCSOLVERTRIANGULAR_CONFIG_trsm_splitfactor_gemm_spdn_param");
    set_by_env(cfg.herk_tri_cfg.partition_algorithm, "ESPRESO_SCSOLVERTRIANGULAR_CONFIG_trsm_splitfactor_gemm_spdn_param");
    set_by_env(cfg.herk_tri_cfg.partition_parameter, "ESPRESO_SCSOLVERTRIANGULAR_CONFIG_trsm_splitfactor_gemm_spdn_param");
}



template<typename T, typename I>
SchurComplementSolver<T,I>::SchurComplementSolver()
{
    ext = std::make_unique<Schur_Complement_Solver_External_Representation<T,I>>();

    init_config<T,I>(ext->cfg);

    ext->A11_solver_factors = DirectSparseSolver<T,I>::factorsSymmetry();
    ext->A11_solver_factors_get_L = ((ext->A11_solver_factors == Solver_Factors::HERMITIAN_LOWER) || (ext->A11_solver_factors == Solver_Factors::NONSYMMETRIC_BOTH));
    ext->A11_solver_factors_get_U = ((ext->A11_solver_factors == Solver_Factors::HERMITIAN_UPPER) || (ext->A11_solver_factors == Solver_Factors::NONSYMMETRIC_BOTH));
    ext->system_is_hermitian = ((ext->A11_solver_factors == Solver_Factors::HERMITIAN_LOWER) || (ext->A11_solver_factors == Solver_Factors::HERMITIAN_UPPER));
    if(ext->A11_solver_factors == Solver_Factors::NONSYMMETRIC_BOTH) {
        eslog::error("nonsymmetric is not supported yet\n");
    }
}



template<typename T, typename I>
SchurComplementSolver<T,I>::~SchurComplementSolver()
{
}



template<typename T, typename I>
SchurComplementSolver<T,I>::SchurComplementSolver(SchurComplementSolver && other) = default;



template<typename T, typename I>
SchurComplementSolver<T,I> & SchurComplementSolver<T,I>::operator=(SchurComplementSolver && other) = default;



template <typename T, typename I>
const char * SchurComplementSolver<T,I>::name()
{
    return "TRIANGULAR";
}



template<typename T, typename I>
void SchurComplementSolver<T,I>::commitMatrix(const Matrix_CSR<T,I> & A, I sc_size)
{
    if(ext->stage != 0) eslog::error("SchurComplementSolver::commitMatrix(): wrong stage\n");
    if(A.nrows != A.ncols) eslog::error("matrix has to be square\n");

    ext->A_whole = &A;
    ext->sc_size = sc_size;
    ext->other_size = A.nrows * sc_size;

    extract_submatrices(*ext->A_whole, ext->A11_my, ext->A12_my, ext->A21_my, ext->A22_my, ext->sc_size);

    ext->A11 = &ext->A11_my;
    ext->A12 = &ext->A12_my;
    ext->A21 = &ext->A21_my;
    ext->A22 = &ext->A22_my;
    
    ext->stage = 1;
}



template<typename T, typename I>
void SchurComplementSolver<T,I>::commitMatrix(const Matrix_CSR<T,I> & A11, const Matrix_CSR<T,I> & A12, const Matrix_CSR<T,I> & A21, const Matrix_CSR<T,I> & A22)
{
    if(ext->stage != 0) eslog::error("SchurComplementSolver::commitMatrix(): wrong stage\n");

    ext->A11 = &A11;
    ext->A12 = &A12;
    ext->A21 = &A21;
    ext->A22 = &A22;
    
    ext->stage = 1;
}



template<typename T, typename I>
void SchurComplementSolver<T,I>::factorizeSymbolic()
{
    if(ext->stage != 1) eslog::error("SchurComplementSolver::factorizeSymbolic(): wrong stage\n");

    if(utils::is_complex<T>()) eslog::error("I dont support complex systems\n"); // todo, will need to add conjugations to some places

    ext->A22_new = MatrixCsxView_new<T,I>::from_old(*ext->A22);

    ext->A11_solver.commit(*ext->A11);
    ext->A11_solver.symbolicFactorization();
    ext->A11_solver.getPermutation(ext->perm);
    ext->perm_new = PermutationView_new<I>::from_old(ext->perm);

    if(ext->A11_solver_factors_get_L) {
        ext->A11_solver.getFactorL(ext->L11, true, false);
        ext->L_new = MatrixCsxView_new<T,I>::from_old(ext->L11);
        ext->L_new.prop.uplo = 'L';
        ext->L_new.prop.diag = 'N';
    }
    if(ext->A11_solver_factors_get_U) {
        ext->A11_solver.getFactorU(ext->U11, true, false);
        ext->U_new = MatrixCsxView_new<T,I>::from_old(ext->U11);
        ext->U_new.prop.uplo = 'U';
        ext->U_new.prop.diag = 'N';
    }

    if(ext->system_is_hermitian) {
        if(ext->A11_solver_factors_get_U) {
            ext->L_new = ext->U_new.get_transposed_reordered_view();
            // todo complex: also conjugation
        }

        {
            // determine which of A12 or A21 stores the B matrix based on their nnz
            if(ext->A12->nnz != 0) {
                ext->B_right_input = MatrixCsxView_new<T,I>::from_old(*ext->A12);
            }
            if(ext->A21->nnz != 0) {
                ext->B_right_input = MatrixCsxView_new<T,I>::from_old(*ext->A21).get_transposed_reordered_view();
            }

            MatrixCsxData_new<T,I> Bt_perm;
            Bt_perm.set(ext->B_right_input.nrows, ext->B_right_input.ncols, ext->B_right_input.nnz, ext->B_right_input.order, AllocatorCPU_new::get_singleton());
            Bt_perm.alloc();
            math::operations::permute_csx_csx<T,I>::do_all(&ext->B_right_input, &Bt_perm, &ext->perm_new, nullptr);

            VectorDenseData_new<I> Bt_perm_colpivots;
            Bt_perm_colpivots.set(Bt_perm.ncols, AllocatorCPU_new::get_singleton());
            Bt_perm_colpivots.alloc();
            math::operations::pivots_trails_csx<T,I>::do_all(&Bt_perm, &Bt_perm_colpivots, 'C', 'P', 'N');

            ext->perm_to_sort.set(Bt_perm.ncols, AllocatorCPU_new::get_singleton());
            ext->perm_to_sort.alloc();

            math::operations::sorting_permutation<I,I>::do_all(&Bt_perm_colpivots, &ext->perm_to_sort);

            ext->B_right.set(Bt_perm.nrows, Bt_perm.ncols, Bt_perm.nnz, Bt_perm.order, AllocatorCPU_new::get_singleton());
            ext->B_right.alloc();
            
            math::operations::permute_csx_csx<T,I>::do_all(&Bt_perm, &ext->B_right, nullptr, &ext->perm_to_sort);

            Bt_perm_colpivots.clear();
            Bt_perm.clear();
        }

        ext->X.set(ext->B_right.nrows, ext->B_right.ncols, ext->cfg.order_X, AllocatorCPU_new::get_singleton());

        ext->L_new_to_use = &ext->L_new;
        if(ext->L_new.order != ext->cfg.order_L) {
            ext->L_new_reordered.set(ext->L_new.nrows, ext->L_new.ncols, ext->L_new.nnz, ext->cfg.order_L, AllocatorCPU_new::get_singleton());
            ext->L_new_reordered.alloc();
            ext->op_reorder_L.set_matrix_src(&ext->L_new);
            ext->op_reorder_L.set_matrix_src(&ext->L_new_reordered);
            ext->op_reorder_L.perform_pattern();
            ext->L_new_to_use = &ext->L_new_reordered;
        }

        ext->op_trsm_tri_L.set_config(ext->cfg.trsm_tri_cfg);
        ext->op_trsm_tri_L.set_L(ext->L_new_to_use);
        ext->op_trsm_tri_L.set_X(&ext->X);
        ext->op_trsm_tri_L.set_X_pattern(ext->B_right);
        ext->op_trsm_tri_L.preprocess();

        ext->op_herk_tri.set_config(ext->cfg.herk_tri_cfg);
        ext->op_herk_tri.set_matrix_A(&ext->X);
        // ext->op_herk_tri.set_matrix_C(); // do it in perform
        ext->op_herk_tri.set_coefficients(T{1}, T{0});
        ext->op_herk_tri.set_mode(math::blas::herk_mode::AhA);
        ext->op_herk_tri.set_A_pattern(ext->B_right);
        ext->op_herk_tri.preprocess();
    }
    else {
        eslog::error("not supported yet, todo\n");
    }

    ext->stage = 2;
}



template<typename T, typename I>
void SchurComplementSolver<T,I>::updateMatrixValues()
{
    if(ext->stage < 2) eslog::error("SchurComplementSolver::updateMatrixValues(): wrong stage\n");

    if(ext->A_whole != nullptr) {
        extract_submatrices(*ext->A_whole, ext->A11_my, ext->A12_my, ext->A21_my, ext->A22_my, ext->sc_size);
    }

    if(ext->system_is_hermitian) {
        math::operations::permute_csx_csx<T,I>::do_all(&ext->B_right_input, &ext->B_right, &ext->perm_new, &ext->perm_to_sort);
    
        if(ext->L_new.order != ext->cfg.order_L) {
            ext->op_reorder_L.perform_values();
        }
    }
    else {
        eslog::error("not supported yet, todo\n");
    }
    
    ext->stage = 3;
}



template<typename T, typename I>
template<typename A>
void SchurComplementSolver<T,I>::factorizeNumericAndGetSc(Matrix_Dense<T,I,A> & sc, char uplo, T alpha)
{
    if(ext->stage != 3) eslog::error("SchurComplementSolver::factorizeNumericAndGetSc(): wrong stage\n");
    if(sc.nrows != sc.ncols) eslog::error("sc matrix must be square\n");
    if(sc.nrows != ext->sc_size) eslog::error("sc matrix does not match the previously given sc size\n");

    ext->A11_solver.numericalFactorization();
    if(ext->A11_solver_factors_get_L) {
        ext->A11_solver.getFactorL(ext->L11, false, true);
    }
    if(ext->A11_solver_factors_get_U) {
        ext->A11_solver.getFactorU(ext->U11, false, true);
    }

    if(ext->system_is_hermitian) {
        ext->X.alloc();
        math::operations::convert_csx_dny<T,I>::do_all(&ext->B_right, &ext->X);

        ext->op_trsm_tri_L.perform();

        MatrixDenseData_new<T> sc_tmp1;
        sc_tmp1.set(ext->sc_size, ext->sc_size, 'R', AllocatorCPU_new::get_singleton());
        sc_tmp1.prop.uplo = uplo;
        sc_tmp1.alloc();

        ext->op_herk_tri.set_matrix_C(&sc_tmp1);
        ext->op_herk_tri.perform();

        ext->X.free();

        math::operations::complete_dnx<T>::do_all(&sc_tmp1, true);

        MatrixDenseData_new<T> sc_tmp2;
        sc_tmp2.set(ext->sc_size, ext->sc_size, 'R', AllocatorCPU_new::get_singleton());
        sc_tmp2.prop.uplo = uplo;
        sc_tmp2.alloc();

        math::operations::permute_dnx_dnx<T,I>::do_all(&sc_tmp1, &sc_tmp2, &ext->perm_to_sort, &ext->perm_to_sort);

        sc_tmp1.free();

        MatrixDenseView_new<T> sc_new = MatrixDenseView_new<T>::from_old(sc);
        sc_new.prop.uplo = uplo;

        math::operations::convert_csx_dny<T,I>::do_all(&ext->A22_new, &sc_new);

        math::operations::lincomb_matrix_dnx<T>::do_all(&sc_new, alpha, &sc_new, -alpha, &sc_tmp2);

        sc_tmp2.free();
    }
    else {
        eslog::error("not supported yet, todo\n");
    }

    ext->stage = 4;
}



template<typename T, typename I>
void SchurComplementSolver<T,I>::solveA11(const Vector_Dense<T,I> & rhs, Vector_Dense<T,I> & sol)
{
    if(ext->stage != 4) eslog::error("SchurComplementSolver::factorizeNumericAndGetSc(): wrong stage\n");
    
    ext->A11_solver.solve(rhs, sol);
}



}

#include "math/wrappers/math.sc_solver.inst.hpp"

#endif
