
#include "math/operations/gemm_csx_dny_dny_prune.h"

#include "math/operations/submatrix_dnx_dnx_noncontig.h"
#include "math/operations/supermatrix_dnx_dnx_noncontig.h"
#include "basis/utilities/stacktimer.h"



namespace espreso {
namespace math {
namespace operations {



template<typename T, typename I>
gemm_csx_dny_dny_prune<T,I>::~gemm_csx_dny_dny_prune()
{
    finalize();
}



template<typename T, typename I>
void gemm_csx_dny_dny_prune<T,I>::set_config(char spdn_A_, bool prune_rows_, bool prune_cols_)
{
    spdn_A = spdn_A_;
    prune_rows = prune_rows_;
    prune_cols = prune_cols_;

    set_config_called = true;
}



template<typename T, typename I>
void gemm_csx_dny_dny_prune<T,I>::set_matrix_A(MatrixCsxView_new<T,I> * A_)
{
    A = A_;
}



template<typename T, typename I>
void gemm_csx_dny_dny_prune<T,I>::set_matrix_B(MatrixDenseView_new<T> * B_)
{
    B = B_;
}



template<typename T, typename I>
void gemm_csx_dny_dny_prune<T,I>::set_matrix_C(MatrixDenseView_new<T> * C_)
{
    C = C_;
}



template<typename T, typename I>
void gemm_csx_dny_dny_prune<T,I>::set_coefficients(T alpha_, T beta_)
{
    alpha = alpha_;
    beta = beta_;
}



template<typename T, typename I>
void gemm_csx_dny_dny_prune<T,I>::preprocess()
{
    if(!set_config_called) eslog::error("config is not set\n");
    if(A == nullptr) eslog::error("A is not set\n");
    if(B == nullptr) eslog::error("B is not set\n");
    if(C == nullptr) eslog::error("C is not set\n");
    if(A->nrows != C->nrows || B->ncols != C->ncols || A->ncols != B->nrows) eslog::error("incompatible matrices\n");
    if(B->order != C->order) eslog::error("B and C order must match\n");

    stacktimer::push("gemm_csx_dny_dny_prune::preprocess");

    op_prune_A.set_matrix_src(A);
    op_prune_A.set_pruning_mode(prune_rows, prune_cols);
    op_prune_A.setup();

    m = op_prune_A.get_dst_matrix_nrows();
    n = B->ncols;
    k = op_prune_A.get_dst_matrix_ncols();

    // stacktimer::info("m %zu n %zu k %zu orderA %c orderBC %c spdnA %c", m, n, k, A->order, B->order, spdn_A);

    if(prune_rows) {
        pruned_rows.set(m, AllocatorCPU_new::get_singleton());
        pruned_rows.alloc();
        op_prune_A.set_vector_pruned_rows(&pruned_rows);
    }
    if(prune_cols) {
        pruned_cols.set(k, AllocatorCPU_new::get_singleton());
        pruned_cols.alloc();
        op_prune_A.set_vector_pruned_cols(&pruned_cols);
    }

    op_prune_A.preprocess2();

    if(spdn_A == 'S') {
        A_pruned_sp.set(m, k, A->nnz, A->order, AllocatorCPU_new::get_singleton());

        op_prune_A.set_matrix_dst_sp(&A_pruned_sp);
    }
    if(spdn_A == 'D') {
        A_pruned_dn.set(m, k, A->order, AllocatorCPU_new::get_singleton());

        op_prune_A.set_matrix_dst_dn(&A_pruned_dn);
    }
    B_pruned.set(k, n, B->order, AllocatorCPU_new::get_singleton());
    C_pruned.set(m, n, C->order, AllocatorCPU_new::get_singleton());

    B_to_use = (prune_cols ? &B_pruned : B);
    C_to_use = (prune_rows ? &C_pruned : C);
    
    if(spdn_A == 'S') {
        op_gemm_sp.set_matrix_A(&A_pruned_sp);
        op_gemm_sp.set_matrix_B(B_to_use);
        op_gemm_sp.set_matrix_C(C_to_use);
        op_gemm_sp.set_coefficients(alpha, beta);
    }
    if(spdn_A == 'D') {
        op_gemm_dn.set_matrix_A(&A_pruned_dn);
        op_gemm_dn.set_matrix_B(B_to_use);
        op_gemm_dn.set_matrix_C(C_to_use);
        op_gemm_dn.set_coefficients(alpha, beta);
    }

    stacktimer::pop();

    preprocess_called = true;
}



template<typename T, typename I>
void gemm_csx_dny_dny_prune<T,I>::perform()
{
    if(!preprocess_called) eslog::error("preprocess was not called\n");

    stacktimer::push("gemm_csx_dny_dny_prune::perform");
    stacktimer::info("m %zu n %zu k %zu orderA %c orderBC %c spdnA %c", m, n, k, A->order, B->order, spdn_A);

    if(spdn_A == 'S') {
        A_pruned_sp.alloc();
    }
    if(spdn_A == 'D') {
        A_pruned_dn.alloc();
    }
    op_prune_A.perform();

    if(B_to_use == &B_pruned) {
        B_pruned.alloc();
        submatrix_dnx_dnx_noncontig<T,I>::do_all(B, &B_pruned, &pruned_cols, nullptr);
    }
    if(C_to_use == &C_pruned) {
        C_pruned.alloc();
        submatrix_dnx_dnx_noncontig<T,I>::do_all(C, &C_pruned, &pruned_rows, nullptr);
    }

    if(spdn_A == 'S') {
        op_gemm_sp.perform();
    }
    if(spdn_A == 'D') {
        op_gemm_dn.perform();
    }

    if(B_to_use == &B_pruned) {
        B_pruned.free();
    }
    if(C_to_use == &C_pruned) {
        supermatrix_dnx_dnx_noncontig<T,I>::do_all(&C_pruned, C, &pruned_rows, nullptr);
        C_pruned.free();
    }
    
    if(spdn_A == 'S') {
        A_pruned_sp.free();
    }
    if(spdn_A == 'D') {
        A_pruned_dn.free();
    }

    stacktimer::pop();
}



template<typename T, typename I>
void gemm_csx_dny_dny_prune<T,I>::finalize()
{
    if(preprocess_called) {
        pruned_rows.clear();
        pruned_cols.clear();
        if(spdn_A == 'S') {
            A_pruned_sp.clear();
        }
    }
    preprocess_called = false;
}



#define INSTANTIATE_T_I(T,I) \
template class gemm_csx_dny_dny_prune<T,I>;

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
