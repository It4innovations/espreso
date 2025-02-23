
#include "math/operations/trsm_csx_dny.h"



template<typename T, typename I>
trsm_csx_dny<T,I>::~trsm_csx_dny()
{
    if(preprocess_called)
}



template<typename T, typename I>
void trsm_csx_dny<T,I>::set_system_matrix(MatrixCsxView_new<T,I> * A_)
{
    M = M_;
}



template<typename T, typename I>
void trsm_csx_dny<T,I>::set_rhs_sol(MatrixDenseView_new<T> * X_)
{
    X = X_;
}



template<typename T, typename I>
void trsm_csx_dny<T,I>::preprocess()
{
    if(A == nullptr) eslog::error("matrix A is not set\n");
    if(X == nullptr) eslog::error("matrix X is not set\n");
    if(preprocess_called) eslog::error("preprocess has already been called\n");
    if(A->prop.uplo != 'U' && A->prop.uplo != 'L') eslog::error("invalid A uplo\n");
    if(A->prop.diag != 'U' && A->prop.diag != 'N') eslog::error("invalid A diag\n");

    Y.set(X->nrows, X->ncols, X->order, &AllocatorCPU_new::get_singleton());

    math::spblas::trsm(*A, *X, Y, handle, 'P');
}



template<typename T, typename I>
void trsm_csx_dny<T,I>::perform()
{
    if(A == nullptr) eslog::error("matrix A is not set\n");
    if(X == nullptr) eslog::error("matrix X is not set\n");
    if(!preprocess_called) eslog::error("preprocess has not been called\n");
    if(A->prop.uplo != 'U' && A->prop.uplo != 'L') eslog::error("invalid A uplo\n");
    if(A->prop.diag != 'U' && A->prop.diag != 'N') eslog::error("invalid A diag\n");

    Y_tmp.alloc();

    math::spblas::trsm(*A, *X, Y, handle, 'C');

    copy_matrix_dense<T>::do_all(X, &Y);

    Y_tmp.free();
}



template<typename T, typename I>
void trsm_csx_dny<T,I>::finalize()
{
    if(A == nullptr) eslog::error("matrix A is not set\n");
    if(X == nullptr) eslog::error("matrix X is not set\n");
    if(!preprocess_called) eslog::error("preprocess has not been called\n");

    math::spblas::trsm(*A, *X, Y, handle, 'F');
}
