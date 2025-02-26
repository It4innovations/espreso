
#include "math/operations/auxiliary/trsm_trirhs_chunk_splitrhs.h"

#include "math/operations/submatrix_csx_dny.h"
#include "math/operations/trsm_dnx_dny.h"



namespace espreso {
namespace math {
namespace operations {



template<typename T, typename I>
trsm_trirhs_chunk_splitrhs<T,I>::~trsm_trirhs_chunk_splitrhs()
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
void trsm_trirhs_chunk_splitrhs<T,I>::set_range(size_t rhs_start_, size_t rhs_end_)
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
void trsm_trirhs_chunk_splitrhs<T,I>::set_X_colpivots(VectorDenseView_new<size_t> * X_colpivots_)
{
    X_colpivots = X_colpivots_;
}



template<typename T, typename I>
void trsm_trirhs_chunk_splitrhs<T,I>::preprocess()
{
    if(preprocess_called) eslog::error("preprocess was already called\n");
    if(!set_config_called) eslog::error("config is not set\n");
    if(!set_range_called) eslog::error("range is not set\n");
    if(L == nullptr) eslog::error("matrix L is not set\n");
    if(X == nullptr) eslog::error("matrix X is not set\n");
    if(X_colpivots == nullptr) eslog::error("B colpivots is not set\n");
    if(L->nrows != L->ncols) eslog::error("matrix L must be square\n");
    if(L->nrows != X->nrows) eslog::error("incompatible matrix sizes\n");
    if(L->prop.uplo != 'L') eslog::error("matrix L must have uplo=L\n");

    k_start = X_colpivots->vals[rhs_start];
    k_end = L->nrows;
    k_size = k_end - k_start;

    sub_X = X->get_submatrix_view(k_start, k_end, rhs_start, rhs_size);

    if(cfg.factor_spdn == 'S') {
        op_submatrix_L_sp.set_matrix_src(L);
        op_submatrix_L_sp.set_matrix_dst(&sub_L_sp);
        op_submatrix_L_sp.set_bounds(k_start, k_end, k_start, k_end);
        op_submatrix_L_sp.setup();
        size_t nnz = op_submatrix_L_sp.get_output_matrix_nnz();

        sub_L_sp.set(k_size, k_size, nnz, cfg.factor_order, AllocatorCPU_new::get_singleton());
        sub_L_sp.prop.diag = L->prop.diag;
        sub_L_sp.prop.uplo = L->prop.uplo;

        op_trsm_sp.set_system_matrix(&sub_L_sp);
        op_trsm_sp.set_rhs_sol(&sub_X);
        
        sub_L_sp.alloc();
        op_submatrix_L_sp.perform();
        op_trsm_sp.preprocess();
        sub_L_sp.free();
    }
    if(cfg.factor_spdn == 'D') {
        sub_L_dn.set(k_size, k_size, cfg.factor_order, AllocatorCPU_new::get_singleton());
        sub_L_dn.prop.diag = L->prop.diag;
        sub_L_dn.prop.uplo = L->prop.uplo;
    }

    preprocess_called = true;
}



template<typename T, typename I>
void trsm_trirhs_chunk_splitrhs<T,I>::perform()
{
    if(!preprocess_called) eslog::error("preprocess was not called\n");

    if(cfg.factor_spdn == 'S') {
        sub_L_sp.alloc();
        op_submatrix_L_sp.perform();
        op_trsm_sp.perform();
        sub_L_sp.free();
    }
    if(cfg.factor_spdn == 'D') {
        sub_L_dn.alloc();
        submatrix_csx_dny<T,I>::do_all(L, &sub_L_dn, k_start, k_end, k_start, k_end);
        trsm_dnx_dny<T>::do_all(&sub_L_dn, &sub_X);
        sub_L_dn.free();
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



#define INSTANTIATE_T_I(T,I) \
template class trsm_trirhs_chunk_splitrhs<T,I>;

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
