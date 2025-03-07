
#include "math/operations/auxiliary/trsm_trirhs_chunk_splitfactor.h"

#include "basis/utilities/stacktimer.h"



namespace espreso {
namespace math {
namespace operations {



template<typename T, typename I>
trsm_trirhs_chunk_splitfactor<T,I>::~trsm_trirhs_chunk_splitfactor()
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
void trsm_trirhs_chunk_splitfactor<T,I>::set_X_rowtrails(VectorDenseView_new<I> * X_rowtrails_)
{
    X_rowtrails = X_rowtrails_;
}



template<typename T, typename I>
void trsm_trirhs_chunk_splitfactor<T,I>::preprocess()
{
    if(preprocess_called) eslog::error("preprocess was already called\n");
    if(!set_config_called) eslog::error("config is not set\n");
    if(!set_range_called) eslog::error("range is not set\n");
    if(L == nullptr) eslog::error("matrix L is not set\n");
    if(X == nullptr) eslog::error("matrix X is not set\n");
    if(X_rowtrails == nullptr) eslog::error("B rowtrails is not set\n");

    stacktimer::push("trsm_trirhs_chunk_splitfactor::preprocess");

    size_t rhs_start = 0;
    size_t rhs_end = X_rowtrails->vals[k_end - 1] + 1;
    size_t rhs_size = rhs_end - rhs_start;

    sub_X_top.set_view(k_size, rhs_size, X->ld, X->order, nullptr);
    op_submatrix_X_top.set_matrix_src(X);
    op_submatrix_X_top.set_matrix_dst(&sub_X_top);
    op_submatrix_X_top.set_bounds(k_start, k_end, rhs_start, rhs_end);

    sub_X_bot.set_view(X->nrows - k_end, rhs_size, X->ld, X->order, nullptr);
    op_submatrix_X_bot.set_matrix_src(X);
    op_submatrix_X_bot.set_matrix_dst(&sub_X_bot);
    op_submatrix_X_bot.set_bounds(k_end, X->nrows, rhs_start, rhs_end);

    if(cfg.trsm_factor_spdn == 'S') {
        op_submatrix_L_top_sp.set_matrix_src(L);
        op_submatrix_L_top_sp.set_matrix_dst(&sub_L_top_sp);
        op_submatrix_L_top_sp.set_bounds(k_start, k_end, k_start, k_end);
        op_submatrix_L_top_sp.setup();
        size_t nnz = op_submatrix_L_top_sp.get_output_matrix_nnz();

        sub_L_top_sp.set(k_size, k_size, nnz, cfg.trsm_factor_order, AllocatorCPU_new::get_singleton());
        sub_L_top_sp.prop.diag = L->prop.diag;
        sub_L_top_sp.prop.uplo = L->prop.uplo;

        op_trsm_sp.set_system_matrix(&sub_L_top_sp);
        op_trsm_sp.set_rhs_matrix(&sub_X_top);
        op_trsm_sp.set_solution_matrix(&sub_X_top);
    }
    if(cfg.trsm_factor_spdn == 'D') {
        op_submatrix_L_top_dn.set_matrix_src(L);
        op_submatrix_L_top_dn.set_matrix_dst(&sub_L_top_dn);
        op_submatrix_L_top_dn.set_bounds(k_start, k_end, k_start, k_end);

        sub_L_top_dn.set(k_size, k_size, cfg.trsm_factor_order, AllocatorCPU_new::get_singleton());
        sub_L_top_dn.prop.diag = L->prop.diag;
        sub_L_top_dn.prop.uplo = L->prop.uplo;

        op_trsm_dn.set_system_matrix(&sub_L_top_dn);
        op_trsm_dn.set_rhs_sol(&sub_X_top);
    }

    if(sub_X_bot.nrows > 0) {
        if(cfg.gemm_factor_prune != 'N' || cfg.gemm_factor_spdn == 'S') {
            op_submatrix_L_bot_sp.set_matrix_src(L);
            op_submatrix_L_bot_sp.set_matrix_dst(&sub_L_bot_sp);
            op_submatrix_L_bot_sp.set_bounds(k_end, L->nrows, k_start, k_end);
            op_submatrix_L_bot_sp.setup();
            size_t nnz = op_submatrix_L_bot_sp.get_output_matrix_nnz();

            sub_L_bot_sp.set(sub_X_bot.nrows, k_size, nnz, cfg.gemm_factor_order, AllocatorCPU_new::get_singleton());
        }

        if(cfg.gemm_factor_prune == 'N' && cfg.gemm_factor_spdn == 'S') {
            op_gemm_normal_sp.set_matrix_A(&sub_L_bot_sp);
            op_gemm_normal_sp.set_matrix_B(&sub_X_top);
            op_gemm_normal_sp.set_matrix_C(&sub_X_bot);
            op_gemm_normal_sp.set_coefficients(T{-1}, T{1});
        }
        if(cfg.gemm_factor_prune == 'N' && cfg.gemm_factor_spdn == 'D') {
            op_submatrix_L_bot_dn.set_matrix_src(L);
            op_submatrix_L_bot_dn.set_matrix_dst(&sub_L_bot_dn);
            op_submatrix_L_bot_dn.set_bounds(k_end, L->nrows, k_start, k_end);

            sub_L_bot_dn.set(sub_X_bot.nrows, k_size, cfg.gemm_factor_order, AllocatorCPU_new::get_singleton());

            op_gemm_normal_dn.set_matrix_A(&sub_L_bot_dn);
            op_gemm_normal_dn.set_matrix_B(&sub_X_top);
            op_gemm_normal_dn.set_matrix_C(&sub_X_bot);
            op_gemm_normal_dn.set_coefficients(T{-1}, T{1});
            op_gemm_normal_dn.set_conj(false, false);
        }
        if(cfg.gemm_factor_prune != 'N') {
            bool prune_rows = (cfg.gemm_factor_prune == 'R') || (cfg.gemm_factor_prune == 'A');
            bool prune_cols = (cfg.gemm_factor_prune == 'C') || (cfg.gemm_factor_prune == 'A');
            op_gemm_prune.set_config(cfg.gemm_factor_spdn, prune_rows, prune_cols);
            op_gemm_prune.set_matrix_A(&sub_L_bot_sp);
            op_gemm_prune.set_matrix_B(&sub_X_top);
            op_gemm_prune.set_matrix_C(&sub_X_bot);
            op_gemm_prune.set_coefficients(T{-1}, T{1});

            sub_L_bot_sp.alloc();
            op_submatrix_L_bot_sp.perform();
            op_gemm_prune.preprocess();
            sub_L_bot_sp.free();
        }
    }

    stacktimer::pop();

    preprocess_called = true;
}



template<typename T, typename I>
void trsm_trirhs_chunk_splitfactor<T,I>::perform()
{
    if(!preprocess_called) eslog::error("preprocess was not called\n");

    stacktimer::push("trsm_trirhs_chunk_splitfactor::perform");

    op_submatrix_X_top.perform();
    op_submatrix_X_bot.perform();

    if(cfg.trsm_factor_spdn == 'S') {
        sub_L_top_sp.alloc();
        op_submatrix_L_top_sp.perform();
        op_trsm_sp.perform();
        sub_L_top_sp.free();
    }
    if(cfg.trsm_factor_spdn == 'D') {
        sub_L_top_dn.alloc();
        op_submatrix_L_top_dn.perform_all();
        op_trsm_dn.perform();
        sub_L_top_dn.free();
    }

    if(sub_X_bot.nrows > 0) {
        if(cfg.gemm_factor_prune == 'N' && cfg.gemm_factor_spdn == 'S') {
            sub_L_bot_sp.alloc();
            op_submatrix_L_bot_sp.perform();
            op_gemm_normal_sp.perform();
            sub_L_bot_sp.free();
        }
        if(cfg.gemm_factor_prune == 'N' && cfg.gemm_factor_spdn == 'D') {
            sub_L_bot_dn.alloc();
            op_submatrix_L_bot_dn.perform_all();
            op_gemm_normal_dn.perform();
            sub_L_bot_dn.free();
        }
        if(cfg.gemm_factor_prune != 'N') {
            sub_L_bot_sp.alloc();
            op_submatrix_L_bot_sp.perform();
            op_gemm_prune.perform();
            sub_L_bot_sp.free();
        }
    }

    stacktimer::pop();
}



template<typename T, typename I>
void trsm_trirhs_chunk_splitfactor<T,I>::finalize()
{
    if(preprocess_called) {
        if(sub_X_bot.nrows > 0) {
            if(cfg.gemm_factor_prune != 'N') {
                op_gemm_prune.finalize();
            }
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
