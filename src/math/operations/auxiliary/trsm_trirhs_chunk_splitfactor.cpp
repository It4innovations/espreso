
#include "math/operations/auxiliary/trsm_trirhs_chunk_splitfactor.h"



namespace espreso {
namespace math {
namespace operations {



template<typename T, typename I>
void trsm_trirhs_chunk_splitfactor<T,I>::~trsm_trirhs_chunk_splitfactor()
{
    finalize();
}



template<typename T, typename I>
void trsm_trirhs_chunk_splitfactor<T,I>::set_config(config cfg_)
{
    cfg = cfg_;

    set_config_called = true;
}



template<typename T, typename I>
void trsm_trirhs_chunk_splitfactor<T,I>::set_range(size_t k_start_, size_t k_end_)
{
    k_start = k_start_;
    k_end = k_end_;
    k_size = k_end - k_start;

    set_range_called = true;
}



template<typename T, typename I>
void trsm_trirhs_chunk_splitfactor<T,I>::set_L(MatrixCsxView_new<T,I> * L_)
{
    L = L_;
}



template<typename T, typename I>
void trsm_trirhs_chunk_splitfactor<T,I>::set_X(MatrixDenseView_new<T> * X_)
{
    X = X_;
}



template<typename T, typename I>
void trsm_trirhs_chunk_splitfactor<T,I>::set_B_rowtrails(VectorDenseView_new<T> * B_rowtrails_)
{
    B_rowtrails = B_rowtrails_;
}



template<typename T, typename I>
void trsm_trirhs_chunk_splitfactor<T,I>::preprocess()
{
    if(preprocess_called) eslog::error("preprocess was already called\n");
    if(!set_config_called) eslog::error("config is not set\n");
    if(!set_range_called) eslog::error("range is not set\n");
    if(L == nullptr) eslog::error("matrix L is not set\n");
    if(X == nullptr) eslog::error("matrix X is not set\n");
    if(B_rowtrails == nullptr) eslog::error("B rowtrails is not set\n");

    size_t rhs_start = 0;
    size_t rhs_end = B_rowtrails[k_end - 1] + 1;
    size_t rhs_size = rhs_end - rhs_start;
    sub_X_top.set_view(k_size, rhs_size, X->ld, X->order, X->vals + k_start * X->get_stride_row() + rhs_start * X->get_stride_col());
    sub_X_bot.set_view(X->nrows - k_end, rhs_size, X->ld, X->order, X->vals + k_end * X->get_stride_row() + rhs_start * X->get_stride_col());

    if(cfg.trsm_factor_spdn == 'S') {
        op_submatrix_L_top_sp.set_matrix_src(L);
        op_submatrix_L_top_sp.set_matrix_dst(&sub_L_top.sp);
        op_submatrix_L_top_sp.set_bounds(k_start, k_end, k_start, k_end);
        op_submatrix_L_top_sp.setup();
        size_t nnz = op_submatrix_L_top_sp.get_output_matrix_nnz();

        sub_L_top.sp.set(k_size, k_size, nnz, cfg.trsm_factor_order, AllocatorCPU_new::get_singleton());
        sub_L_top.sp.diag = L->diag;
        sub_L_top.sp.uplo = L->uplo;

        op_trsm_sp.set_system_matrix(sub_L_top.sp);
        op_trsm_sp.set_rhs_sol(sub_X_top);

        sub_L_top.sp.alloc();
        op_submatrix_L_top_sp.perform();
        op_trsm_sp.preprocess();
        sub_L_top.sp.free();
    }
    if(cfg.trsm_factor_spdn == 'D') {
        sub_L_top_dn.set(k_size, k_size, cfg.trsm_factor_order, AllocatorCPU_new::get_singleton());
        sub_L_top.dn.diag = L->diag;
        sub_L_top.dn.uplo = L->uplo;
    }

    if(cfg.gemm_factor_prune == 'Y' || cfg.gemm_factor_spdn == 'S') {
        op_submatrix_L_bot_sp.set_matrix_src(L);
        op_submatrix_L_bot_sp.set_matrix_dst(&sub_L_bot.sp);
        op_submatrix_L_bot_sp.set_bounds(k_end, L.nrows, k_start, k_end);
        op_submatrix_L_bot_sp.setup();
        size_t nnz = op_submatrix_L_bot_sp.get_output_matrix_nnz();

        sub_L_bot.sp.set(sub_X_bot.nrows, k_size, nnz, cfg.gemm_factor_order, AllocatorCPU_new::get_singleton());
    }

    if(cfg.gemm_factor_prune == 'N' && cfg.gemm_factor_spdn == 'S') {
        op_gemm.normal_sp.set_matrix_A(sub_L_bot.sp);
        op_gemm.normal_sp.set_matrix_B(sub_X_top);
        op_gemm.normal_sp.set_matrix_C(sub_X_bot);
        op_gemm.normal_sp.set_coefficients(T{-1}, T{1});

        sub_L_bot.alloc();
        op_submatrix_L_bot_sp.perform();
        op_gemm.normal_sp.preprocess();
        sub_L_bot.free();
    }
    if(cfg.gemm_factor_prune == 'N' && cfg.gemm_factor_spdn == 'D') {
        sub_L_bot.set(sub_X_bot.nrows, k_size, cfg.gemm_factor_order, AllocatorCPU_new::get_singleton());
    }
    if(cfg.gemm_factor_prune == 'Y') {
        op_gemm.prune.set_config(cfg.gemm_factor_spdn);
        op_gemm.prune.set_matrix_A(sub_L_bot.sp);
        op_gemm.prune.set_matrix_B(sub_X_top);
        op_gemm.prune.set_matrix_C(sub_X_bot);
        op_gemm.prune.set_coefficients(T{-1}, T{1});
        
        sub_L_bot.alloc();
        op_submatrix_L_bot_sp.perform();
        op_gemm.prune.preprocess();
        sub_L_bot.free();
    }

    preprocess_called = true;
}



template<typename T, typename I>
void trsm_trirhs_chunk_splitfactor<T,I>::perform()
{
    if(!preprocess_called) eslog::error("preprocess was not called\n");

    if(cfg.trsm_factor_spdn == 'S') {
        sub_L_top.sp.alloc();
        op_submatrix_L_top_sp.perform();
        op_trsm_sp.perform();
        sub_L_top.sp.free();
    }
    if(cfg.trsm_factor_spdn == 'D') {
        sub_L_top_dn.alloc();
        submatrix_csx_dny<T,I>::do_all(L, &sub_L_top.dn, k_start, k_end, k_start, k_end);
        trsm_dnx_dny<T,I>::do_all(&sub_L_top.dn, sub_X_top);
        sub_L_top_dn.free();
    }

    if(gemm_factor_prune == 'N' && gemm_factor_spdn == 'S') {
        sub_L_bot.sp.alloc();
        op_submatrix_L_bot_sp.perform();
        op_gemm.normal_sp.perform();
        sub_L_bot.sp.free();
    }
    if(gemm_factor_prune == 'N' && gemm_factor_spdn == 'D') {
        sub_L_bot.dn.alloc();
        submatrix_csx_dny<T,I>::do_all(L, &sub_L_bot.dn, k_end, L.nrows, k_start, k_end);
        gemm_dnx_dny_dnz<T,I>::do_all(&sub_L_bot.dn, sub_X_top, sub_X_bot, T{-1}, T{1});
        sub_L_bot.dn.free();
    }
    if(gemm_factor_prune == 'Y') {
        sub_L_bot.sp.alloc();
        op_submatrix_L_bot_sp.perform();
        op_gemm.prune_sp.perform();
        sub_L_bot.sp.free();
    }
}



template<typename T, typename I>
void trsm_trirhs_chunk_splitfactor<T,I>::finalize()
{
    if(preprocess_called) {
        op_trsm_sp.finalize();
        if(gemm_factor_prune == 'N' && gemm_factor_spdn == 'S') {
            op_gemm.normal_sp.finalize();
        }
        if(gemm_factor_prune == 'Y') {
            op_gemm.prune.finalize();
        }
    }
    preprocess_called = false;
}



#define INSTANTIATE_T_I(T,I) \
template class trsm_trirhs_chunk_splitfactor<T,I>;

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
