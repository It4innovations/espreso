
#include "math/operations/gemm_csx_dny_dny.h"



template<typename T, typename I>
gemm_csx_dny_dny<T,I>::~gemm_csx_dny_dny()
{
    finalize();
}



template<typename T, typename I>
void gemm_csx_dny_dny<T,I>::set_matrix_A(MatrixCsxView_new<T,I> * A_)
{
    A = A_;
}



template<typename T, typename I>
void gemm_csx_dny_dny<T,I>::set_matrix_B(MatrixDenseView_new<T> * B_)
{
    B = B_;
}



template<typename T, typename I>
void gemm_csx_dny_dny<T,I>::set_matrix_C(MatrixDenseView_new<T> * C_)
{
    C = C_;
}



template<typename T, typename I>
void gemm_csx_dny_dny<T,I>::set_coefficients(T alpha_, T beta_)
{
    alpha = alpha_;
    beta = beta_;
}



template<typename T, typename I>
void gemm_csx_dny_dny<T,I>::preprocess()
{
    if(preprocess_called) eslog::error("preproces has already been called\n");
    if(A == nullptr) eslog::error("matrix A is not set\n");
    if(B == nullptr) eslog::error("matrix B is not set\n");
    if(C == nullptr) eslog::error("matrix C is not set\n");
    if(A->nrows != C.nrows || B.ncols != C.ncols || A.ncols != B.nrows) eslog::error("incompatible matrix sizes\n");
    if(B.order != C.order) eslog::error("order of B and C must match\n");

    math::blas::mm(*A, *B, *C, alpha, beta, handle_abc, 'P');

    preprocess_called = true;
}



template<typename T, typename I>
void gemm_csx_dny_dny<T,I>::perform()
{
    if(!preprocess_called) eslog::error("preprocess has not been called\n");
    if(A == nullptr) eslog::error("matrix A is not set\n");
    if(B == nullptr) eslog::error("matrix B is not set\n");
    if(C == nullptr) eslog::error("matrix C is not set\n");
    if(A->nrows != C.nrows || B.ncols != C.ncols || A.ncols != B.nrows) eslog::error("incompatible matrix sizes\n");
    if(B.order != C.order) eslog::error("order of B and C must match\n");

    math::blas::mm(*A, *B, *C, alpha, beta, handle_abc, 'C');
}



template<typename T, typename I>
void gemm_csx_dny_dny<T,I>::finalize()
{
    if(preprocess_called) {
        math::blas::mm(*A, *B, *C, alpha, beta, handle_abc, 'F');
    }
    preprocess_called = false;
}



template<typename T, typename I>
void gemm_csx_dny_dny<T,I>::do_all(MatrixCsxView_new<T,I> * A, MatrixDenseView_new<T> * B, MatrixDenseView_new<T> * C, T alpha, T beta)
{
    gemm_csx_dny_dny<T,I> instance;
    instance.set_matrix_A(A);
    instance.set_matrix_B(B);
    instance.set_matrix_C(C);
    instance.set_coefficients(alpha, beta);
    instance.preprocess();
    instance.perform();
    instance.finalize();
}
