
#include "math/operations/gemm_dnx_dny_dnz.h"

#include "math/wrappers/math.blas.h"



template<typename T>
void gemm_dnx_dny_dnz<T>::set_matrix_A(MatrixDenseView_new<T,I> * A_)
{
    A = A_;
}



template<typename T>
void gemm_dnx_dny_dnz<T>::set_matrix_B(MatrixDenseView_new<T,I> * B_)
{
    B = B_;
}



template<typename T>
void gemm_dnx_dny_dnz<T>::set_matrix_C(MatrixDenseView_new<T,I> * C_)
{
    C = C_;
}



template<typename T>
void gemm_dnx_dny_dnz<T>::set_coefficients(T alpha_, T beta_)
{
    alpha = alpha_;
    beta = beta_;
}



template<typename T>
void gemm_dnx_dny_dnz<T>::perform()
{
    if(A == nullptr) eslog::error("matrix A is not set\n");
    if(B == nullptr) eslog::error("matrix B is not set\n");
    if(C == nullptr) eslog::error("matrix C is not set\n");
    if(A->nrows != C->nrows || B->ncols != C.ncols || A->ncols != B.nrows) eslog::error("incompatible matrices");

    math::blas::gemm(*A, *B, *C, alpha, beta);
}



template<typename T>
void gemm_dnx_dny_dnz<T>::do_all(MatrixDenseView_new<T,I> * A, MatrixDenseView_new<T,I> * B, MatrixDenseView_new<T,I> * C, T alpha, T beta)
{
    gemm_dnx_dny_dnz<T> instance;
    instance.set_matrix_A(A);
    instance.set_matrix_B(B);
    instance.set_matrix_C(C);
    instance.set_coefficients(alpha, beta);
    instance.perform();
}
