
#include "math/operations/auxiliary/trsm_trirhs_chunk_splitrhs.h"



template<typename T, typename I>
void trsm_trirhs_chunk_splitrhs<T,I>::~trsm_trirhs_chunk_splitrhs()
{
    finalize();
}



template<typename T, typename I>
void trsm_trirhs_chunk_splitrhs<T,I>::set_config(config cfg_)
{
    cfg = cfg_;

    set_config_called = true;
}



template<typename T, typename I>
void trsm_trirhs_chunk_splitrhs<T,I>::set_range(size_t rhs_start, size_t rhs_end)
{
    rhs_start = rhs_start_;
    rhs_end = rhs_end_;
    rhs_size = rhs_end - rhs_start;

    set_range_called = true;
}



template<typename T, typename I>
void trsm_trirhs_chunk_splitrhs<T,I>::set_L(MatrixCsxView_new<T,I> * L_)
{
    L = L_;
}



template<typename T, typename I>
void trsm_trirhs_chunk_splitrhs<T,I>::set_X(MatrixDenseView_new<T> * X_)
{
    X = X_;
}



template<typename T, typename I>
void trsm_trirhs_chunk_splitrhs<T,I>::set_B_colpivots(VectorDenseView_new<T> * B_colpivots_)
{
    B_colpivots = B_colpivots_;
}



template<typename T, typename I>
void trsm_trirhs_chunk_splitrhs<T,I>::preprocess()
{
    if(preprocess_called) eslog::error("preprocess was already called\n");
    if(!set_config_called) eslog::error("config is not set\n");
    if(!set_range_called) eslog::error("range is not set\n");
    if(L == nullptr) eslog::error("matrix L is not set\n");
    if(X == nullptr) eslog::error("matrix X is not set\n");
    if(B_colpivots == nullptr) eslog::error("B colpivots is not set\n");
    if(L->nrows != L->ncols) eslog::error("matrix L must be square\n");
    if(L->nrows != B->nrows) eslog::error("incompatible matrix sizes\n");
    if(L->prop.uplo != 'L') eslog::error("matrix L must have uplo=L\n");

    k_start = B_colpivots[rhs_start];
    k_end = L->nrows;
    k_size = k_end - k_start;

    sub_X.set_view(k_size, rhs_size, X->size, X->order, X->vals + k_start * X->get_stride_row() + rhs_start * X->get_stride_col());

    if(cfg.factor_spdn == 'S') {
        op_submatrix_L_sp.set_matrix_src(L);
        op_submatrix_L_sp.set_matrix_dst(&sub_L.sp);
        op_submatrix_L_sp.set_bounds(k_start, k_end, k_start, k_end);
        op_submatrix_L_sp.setup();
        size_t nnz = op_submatrix_L_sp.get_output_matrix_nnz();

        sub_L.sp.set(k_size, k_size, nnz, cfg.factor_order, AllocatorCPU_new::get_singleton());
        sub_L.sp.diag = L->diag;
        sub_L.sp.uplo = L->uplo;

        op_trsm_sp.set_system_matrix(&sub_L.sp);
        op_trsm_sp.set_rhs_sol(&sub_X);
        
        sub_L.sp.alloc();
        op_submatrix_L_sp.perform();
        op_trsm_sp.preprocess();
        sub_L.sp.free();
    }
    if(cfg.factor_spdn == 'D') {
        sub_L.dn.set(k_size, k_size, cfg.factor_order, AllocatorCPU_new::get_singleton());
        sub_L.dn.diag = L->diag;
        sub_L.dn.uplo = L->uplo;
    }

    preprocess_called = true;
}



template<typename T, typename I>
void trsm_trirhs_chunk_splitrhs<T,I>::perform()
{
    if(!preprocess_called) eslog::error("preprocess was not called\n");

    if(cfg.factor_spdn == 'S') {
        sub_L.sp.alloc();
        op_submatrix_L_sp.perform();
        op_trsm_sp.perform();
        sub_L.sp.free();
    }
    if(cfg.factor_spdn == 'D') {
        sub_L.dn.alloc();
        submatrix_csx_dny<T,I>::do_all(L, &sub_L.dn, k_start, k_end, k_start, k_end);
        trsm_dnx_dny<T>::do_all(&sub_L.dn, sub_X);
        sub_L.dn.free();
    }
}



template<typename T, typename I>
void trsm_trirhs_chunk_splitrhs<T,I>::finalize()
{
    if(preprocess_called) {
        op_trsm_sp.finalize();
    }
    preprocess_called = false;
}




