
#include "math/operations/gemm_csx_dny_dny_prune.h"



namespace espreso {
namespace math {
namespace operations {



template<typename T, typename I>
gemm_csx_dny_dny_prune<T,I>::~gemm_csx_dny_dny_prune()
{
    finalize();
}



template<typename T, typename I>
gemm_csx_dny_dny_prune<T,I>::set_config(char spdn_A_)
{
    spdn_A = spdn_A_;

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

    size_t num_nonempty = 0;
    size_t A_size_primary = A.get_primary_size();
    I * A_ptrs = A->ptrs;
    for(size_t ip = 0; ip < A_size_primary; ip++) {
        I start = A_ptrs[ip];
        I end = A_ptrs[ip+1];
        if(start != end) {
            num_nonempty++;
        }
    }

    nonempty.set(num_nonempty, AllocatorCPU_new::get_singleton());

    if(spdn_A == 'S') {
        A_pruned_ptrs.set(num_nonempty + 1, AllocatorCPU_new::get_singleton());
    
        size_t idx_nonempty = 0;
        for(size_t ip = 0; ip < A_size_primary; ip++) {
            I start = A_ptrs[ip];
            I end = A_ptrs[ip+1];
            if(start != end) {
                nonempty.vals[idx_nonempty] = ip;
                A_pruned_ptrs.vals[idx_nonempty] = start;
                idx_nonempty++;
            }
        }
        A_pruned_ptrs.vals[num_nonempty] = A->nnz;
    }
    if(spdn_A == 'D') {
        size_t idx_nonempty = 0;
        for(size_t ip = 0; ip < A_size_primary; ip++) {
            I start = A_ptrs[ip];
            I end = A_ptrs[ip+1];
            if(start != end) {
                nonempty.vals[idx_nonempty] = ip;
                idx_nonempty++;
            }
        }
    }

    size_t m = ((A->order == 'R') ? num_nonempty : A->nrows);
    size_t n = B->ncols;
    size_t k = ((A->order == 'C') ? num_nonempty : A->ncols);

    if(spdn_A == 'S') {
        A_pruned.sp.set_view(m, k, A->nnz, A->order, A_pruned_ptrs.vals, A->idxs, A->vals);
    }
    if(spdn_A == 'D') {
        A_pruned.set(m, k, A->order, AllocatorCPU_new::get_singleton());
    }
    B_pruned.set(k, n, B->order, AllocatorCPU_new::get_singleton());
    C_pruned.set(m, n, C->order, AllocatorCPU_new::get_singleton());

    B_to_use = ((A->order == 'C') ? &B_pruned : B);
    C_to_use = ((A->order == 'R') ? C : &C_pruned);
    
    if(spdn_A == 'S') {
        op_gemm_sp.set_matrix_A(&A_pruned);
        op_gemm_sp.set_matrix_B(B_to_use);
        op_gemm_sp.set_matrix_C(C_to_use);
        op_gemm_sp.set_coefficients(alpha, beta);
        op_gemm_sp.preprocess();
    }

    preprocess_called = true;
}



template<typename T, typename I>
void gemm_csx_dny_dny_prune<T,I>::perform()
{
    if(!preprocess_called) eslog::error("preprocess was not called\n");

    if(spdn_A == 'D') {
        A_pruned.dn.alloc();
        std::fill_n(A_pruned.dn.vals, A_pruned.dn.get_num_blocks() * A_pruned.dn.ld, 0);
        for(size_t idx_nonempty = 0; idx_nonempty < nonempty.size; idx_nonempty++) {
            size_t ip = nonempty[idx_nonempty];
            I start = A_ptrs[ip];
            I end = A_ptrs[ip+1];
            for(I i = start; i < end; i++) {
                I is = A_idxs[i];
                T v = A_vals[i];
                A_pruned.dn.vals[ip * A_pruned.dn.ld + is] = v;
            }
        }
    }

    if(B_to_use == &B_pruned) {
        B_pruned.alloc();
        submatrix_dnx_dnx_noncontig<T>::do_all(B, &B_pruned, &nonempty, nullptr);
    }
    if(C_to_use == &C_pruned) {
        C_pruned.alloc();
        submatrix_dnx_dnx_noncontig<T>::do_all(C, &C_pruned, &nonempty, nullptr);
    }

    if(spdn_A == 'S') {
        op_gemm_sp.perform();
    }
    if(spdn_A == 'D') {
        gemm_dnx_dny_dnz<T>::do_all(&A_pruned.dn, B_to_use, C_to_use, alpha, beta);
    }

    if(B_to_use == &B_pruned) {
        B_pruned.free();
    }
    if(C_to_use == &C_pruned) {
        supermatrix_dnx_dnx_noncontig<T>::do_all(&C_pruned, C, &nonempty, nullptr);
        C_pruned.free();
    }
    
    if(spdn_A == 'D') {
        A_pruned.dn.free();
    }
}



template<typename T, typename I>
void gemm_csx_dny_dny_prune<T,I>::finalize()
{
    if(preprocess_called) {
        if(spdn_A == 'S') {
            op_gemm_sp.finalize();
            A_pruned_ptrs.clear();
        }
        nonempty.clear();
    }
    preprocess_called = false;
}



template<typename T, typename I>
void gemm_csx_dny_dny_prune<T,I>::do_all(MatrixCsxView_new<T,I> * A, MatrixDenseView_new<T> * B, MatrixDenseView_new<T> * C, T alpha, T beta, char spdn_A)
{
    gemm_csx_dny_dny_prune<T,I> instance;
    instance.set_config(spdn_A);
    instance.set_matrix_A(A);
    instance.set_matrix_B(B);
    instance.set_matrix_C(C);
    instance.set_coefficients(alpha, beta);
    instance.preprocess();
    instance.perform();
    instance.finalize();
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
        /* INSTANTIATE_T(std::complex<double>) */

            INSTANTIATE

        #undef INSTANTIATE
    #undef INSTANTIATE_T
#undef INSTANTIATE_T_I



}
}
}
