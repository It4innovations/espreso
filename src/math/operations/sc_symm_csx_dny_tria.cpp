
#include "math/operations/sc_symm_csx_dny_tria.h"

#include "math/primitives_new/allocator_new.h"
#include "math/operations/permute_csx_csx.h"
#include "math/operations/pivots_trails_csx.h"
#include "math/operations/sorting_permutation.h"
#include "math/operations/convert_csx_dny.h"
#include "math/operations/complete_dnx.h"
#include "math/operations/permute_dnx_dnx.h"
#include "math/operations/copy_dnx.h"
#include "basis/utilities/stacktimer.h"



namespace espreso {
namespace math {
namespace operations {



template<typename T, typename I>
void sc_symm_csx_dny_tria<T,I>::set_config(config cfg_)
{
    if(called_set_config) eslog::error("config is already set\n");

    cfg = cfg_;

    called_set_config = true;
}



template<typename T, typename I>
void sc_symm_csx_dny_tria<T,I>::set_coefficients(Treal alpha_)
{
    alpha = alpha_;
}



template<typename T, typename I>
void sc_symm_csx_dny_tria<T,I>::set_A11_solver(DirectSparseSolver<T,I> * A11_solver_)
{
    if(A11_solver != nullptr) eslog::error("A11_solver is already set\n");

    A11_solver = A11_solver_;
}



template<typename T, typename I>
void sc_symm_csx_dny_tria<T,I>::set_A12(MatrixCsxView_new<T,I> * A12_)
{
    if(A12 != nullptr) eslog::error("matrix A12 is already set\n");

    A12 = A12_;
}



template<typename T, typename I>
void sc_symm_csx_dny_tria<T,I>::set_sc(MatrixDenseView_new<T> * sc_)
{
    if(sc != nullptr) eslog::error("matrix sc is already set\n");

    sc = sc_;
}



template<typename T, typename I>
void sc_symm_csx_dny_tria<T,I>::preprocess()
{
    stacktimer::push("sc_symm_csx_dny_tria::preprocess");

    if(!called_set_config) eslog::error("config is not set\n");
    if(A11_solver == nullptr) eslog::error("A11 solver is not set\n");
    if(A12 == nullptr) eslog::error("matrix A12 is not set\n");
    if(sc == nullptr) eslog::error("sc matrix is not set\n");
    if(!A12->ator->is_data_accessible_cpu()) eslog::error("matrix A12 must be cpu-accessible\n");
    if(!sc->ator->is_data_accessible_cpu()) eslog::error("matrix sc must be cpu-accessible\n");
    if((size_t)A11_solver->getMatrixSize() != A12->nrows) eslog::error("incompatible matrices\n");
    if(sc->ncols != A12->ncols) eslog::error("incompatible matrices\n");
    if(sc->prop.uplo != 'L' && sc->prop.uplo != 'U') eslog::error("wrong sc uplo\n");

    solver_factor_uplo = '_';
    if(DirectSparseSolver<T,I>::factorsSymmetry() == Solver_Factors::HERMITIAN_LOWER) solver_factor_uplo = 'L';
    if(DirectSparseSolver<T,I>::factorsSymmetry() == Solver_Factors::HERMITIAN_UPPER) solver_factor_uplo = 'U';
    if(solver_factor_uplo == '_') eslog::error("wrong sparse solver, must be symmetric\n");
    if(!DirectSparseSolver<T,I>::provideFactors()) eslog::error("wrong sparse solver, must provide factors\n");

    I factor_nnz = A11_solver->getFactorNnz();
    A11size = A11_solver->getMatrixSize();

    factor_U_row.set(A11size, A11size, factor_nnz, 'R', AllocatorCPU_new::get_singleton());
    factor_U_row.prop.uplo = 'U';
    factor_U_row.prop.diag = 'N';

    factor_L_row.set(A11size, A11size, factor_nnz, 'R', AllocatorCPU_new::get_singleton());
    factor_L_row.prop.uplo = 'L';
    factor_L_row.prop.diag = 'N';

    need_reorder_factor_L2U = ((solver_factor_uplo == 'L') && (cfg.order_L == 'C'));
    need_reorder_factor_U2L = ((solver_factor_uplo == 'U') && (cfg.order_L == 'R'));

    if(solver_factor_uplo == 'U' || need_reorder_factor_L2U) {
        factor_U_row.alloc();
    }
    if(solver_factor_uplo == 'L' || need_reorder_factor_U2L) {
        factor_L_row.alloc();
    }

    stacktimer::push("extract_factors");
    if(solver_factor_uplo == 'U') {
        Matrix_CSR<T,I> factor_old_U = MatrixCsxView_new<T,I>::template to_old<cpu_allocator>(factor_U_row);
        A11_solver->getFactorU(factor_old_U, true, false);
    }
    if(solver_factor_uplo == 'L') {
        Matrix_CSR<T,I> factor_old_L = MatrixCsxView_new<T,I>::template to_old<cpu_allocator>(factor_L_row);
        A11_solver->getFactorL(factor_old_L, true, false);
    }
    stacktimer::pop();

    L_row = factor_L_row;
    L_col = factor_U_row.get_transposed_reordered_view();
    U_row = factor_U_row;
    U_col = factor_L_row.get_transposed_reordered_view();

    if(need_reorder_factor_L2U) {
        op_L2U.set_matrix_src(&L_row);
        op_L2U.set_matrix_dst(&L_col);
        op_L2U.perform_pattern();
    }
    if(need_reorder_factor_U2L) {
        op_U2L.set_matrix_src(&U_row);
        op_U2L.set_matrix_dst(&U_col);
        op_U2L.perform_pattern();
    }

    if(cfg.order_L == 'R') L_to_use = &L_row;
    if(cfg.order_L == 'C') L_to_use = &L_col;

    X_sp.set(A12->nrows, A12->ncols, A12->nnz, A12->order, AllocatorCPU_new::get_singleton());
    X_sp.alloc();

    {
        MatrixCsxData_new<T,I> X_sp_tmp;
        X_sp_tmp.set(X_sp.nrows, X_sp.ncols, X_sp.nnz, X_sp.order, AllocatorCPU_new::get_singleton());
        X_sp_tmp.alloc();

        {
            perm_fillreduce.set(A12->nrows, AllocatorCPU_new::get_singleton());
            perm_fillreduce.alloc();
            Permutation<I> perm_fillreduce_old = PermutationView_new<I>::template to_old<cpu_allocator>(perm_fillreduce);
            A11_solver->getPermutation(perm_fillreduce_old);

            math::operations::permute_csx_csx<T,I>::do_all(A12, &X_sp_tmp, &perm_fillreduce, nullptr);
        }

        {
            VectorDenseData_new<I> colpivots;
            colpivots.set(A12->ncols, AllocatorCPU_new::get_singleton());
            colpivots.alloc();

            math::operations::pivots_trails_csx<T,I>::do_all(&X_sp_tmp, &colpivots, 'C', 'P', '_');

            perm_to_sort_A12_cols.set(A12->ncols, AllocatorCPU_new::get_singleton());
            perm_to_sort_A12_cols.alloc();

            perm_to_sort_back_sc = perm_to_sort_A12_cols.get_inverse_view();

            math::operations::sorting_permutation<I,I>::do_all(&colpivots, &perm_to_sort_A12_cols);
        }

        math::operations::permute_csx_csx<T,I>::do_all(&X_sp_tmp, &X_sp, nullptr, &perm_to_sort_A12_cols);
    }

    X_dn.set(X_sp.nrows, X_sp.ncols, cfg.order_X, AllocatorCPU_new::get_singleton());

    sc_tmp1.set(sc->nrows, sc->ncols, sc->order, AllocatorCPU_new::get_singleton());
    sc_tmp1.prop.uplo = sc->prop.uplo;

    sc_tmp2.set(sc->nrows, sc->ncols, sc->order, AllocatorCPU_new::get_singleton());
    sc_tmp2.prop.uplo = sc->prop.uplo;

    op_trsm.set_config(cfg.cfg_trsm);
    op_trsm.set_L(L_to_use);
    op_trsm.set_X(&X_dn);
    op_trsm.calc_X_pattern(X_sp);
    op_trsm.preprocess();

    op_herk.set_config(cfg.cfg_herk);
    op_herk.set_matrix_A(&X_dn);
    op_herk.set_matrix_C(&sc_tmp1);
    op_herk.set_coefficients(-alpha, 0 * alpha);  // beta=0, because I assume A22=0
    op_herk.set_mode(blas::herk_mode::AhA);
    op_herk.calc_A_pattern(X_sp);
    op_herk.preprocess();

    stacktimer::pop();

    called_preprocess = true;
}



template<typename T, typename I>
void sc_symm_csx_dny_tria<T,I>::perform()
{
    stacktimer::push("sc_symm_csx_dny_tria::perform");

    stacktimer::push("extract_factors");
    if(solver_factor_uplo == 'U') {
        Matrix_CSR<T,I,gpu::mgm::Ah> factor_old = MatrixCsxView_new<T,I>::template to_old<gpu::mgm::Ah>(factor_U_row);
        A11_solver->getFactorU(factor_old, false, true);
    }
    if(solver_factor_uplo == 'L') {
        Matrix_CSR<T,I,gpu::mgm::Ah> factor_old = MatrixCsxView_new<T,I>::template to_old<gpu::mgm::Ah>(factor_L_row);
        A11_solver->getFactorL(factor_old, false, true);
    }
    stacktimer::pop();

    if(need_reorder_factor_L2U) {
        op_L2U.perform_values();
    }
    if(need_reorder_factor_U2L) {
        op_U2L.perform_values();
    }

    X_dn.alloc();
    sc_tmp1.alloc();
    sc_tmp2.alloc();

    convert_csx_dny<T,I>::do_all(&X_sp, &X_dn);

    op_trsm.perform();

    op_herk.perform();

    complete_dnx<T>::do_all(&sc_tmp1, sc->prop.uplo, true);

    permute_dnx_dnx<T,I>::do_all(&sc_tmp1, &sc_tmp2, &perm_to_sort_back_sc, &perm_to_sort_back_sc);

    copy_dnx<T>::do_all(&sc_tmp2, sc, false);

    X_dn.free();
    sc_tmp1.free();
    sc_tmp2.free();

    stacktimer::pop();
}



#define INSTANTIATE_T_I(T,I) \
template class sc_symm_csx_dny_tria<T,I>;

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

