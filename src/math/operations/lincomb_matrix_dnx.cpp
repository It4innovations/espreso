
#include "math/operations/lincomb_matrix_dnx.h"



namespace espreso {
namespace math {
namespace operations {



template<typename T>
void lincomb_matrix_dnx<T>::set_matrix_X(MatrixDenseView_new<T> * X_)
{
    X = X_;
}



template<typename T>
void lincomb_matrix_dnx<T>::set_matrix_A(MatrixDenseView_new<T> * A_)
{
    A = A_;
}



template<typename T>
void lincomb_matrix_dnx<T>::set_matrix_B(MatrixDenseView_new<T> * B_)
{
    B = B_;
}



template<typename T>
void lincomb_matrix_dnx<T>::set_coefficients(T alpha_, T beta_)
{
    alpha = alpha_;
    beta = beta_;
}



template<typename T>
void lincomb_matrix_dnx<T>::perform()
{
    if(X == nullptr) eslog::error("result matrix X is not set\n");
    if(alpha != 0 && A == nullptr) eslog::error("matrix A is not set\n");
    if(beta != 0 && B == nullptr) eslog::error("matrix B is not set\n");

    if(alpha == 0 && beta == 0) {
        lincomb_matrix_dnx<T>::perform_zero(*X);
    }
    if(alpha != 0 && beta == 0) {
        lincomb_matrix_dnx<T>::perform_one(*X, alpha, *A);
    }
    if(alpha == 0 && beta != 0) {
        lincomb_matrix_dnx<T>::perform_one(*X, beta, *B);
    }
    if(alpha != 0 && beta != 0) {
        lincomb_matrix_dnx<T>::perform_two(*X, alpha, *A, beta, *B);
    }
}



template<typename T>
void lincomb_matrix_dnx<T>::do_all(MatrixDenseView_new<T> * X, T alpha, MatrixDenseView_new<T> * A, T beta, MatrixDenseView_new<T> * B)
{
    lincomb_matrix_dnx<T> instance;
    instance.set_matrix_X(X);
    instance.set_matrix_A(A);
    instance.set_matrix_B(B);
    instance.set_coefficients(alpha, beta);
    instance.perform();
}



template<typename T>
void lincomb_matrix_dnx<T>::perform_zero(MatrixDenseView_new<T> & X)
{
    fill_dnx<T>(&X, T{0});
}



template<typename T>
void lincomb_matrix_dnx<T>::perform_one(MatrixDenseView_new<T> & X, T alpha, MatrixDenseView_new<T> & A)
{
    if(A.nrows != X.nrows || A.ncols != X.ncols) eslog::error("matrix sizes dont match\n");
    if(A.order != X.order) eslog::error("matrix orders dont match\n");

    T * vals_X = X.vals;
    T * vals_A = A.vals;
    size_t size_primary = X.get_size_primary();
    size_t size_secdary = X.get_size_secdary();
    bool move_start = ((X.uplo == 'U' && X.order == 'R') || (X.uplo == 'L' && X.order == 'C'));
    bool move_end   = ((X.uplo == 'L' && X.order == 'R') || (X.uplo == 'U' && X.order == 'C'));
    for(size_t i = 0; i < size_primary; i++) {
        size_t start = 0;
        size_t end = size_secdary;
        if(move_start) start = i;
        if(move_end) end = i;
        size_t size = end - start;
        T * sub_X = vals_X + i * X.ld;
        T * sub_A = vals_A + i * A.ld;
        for(size_t j = start; j < end; j++) {
            sub_X[j] = alpha * sub_A[j];
        }
    }
}



template<typename T>
void lincomb_matrix_dnx<T>::perform_two(MatrixDenseView_new<T> & X, T alpha, MatrixDenseView_new<T> & A, T beta, MatrixDenseView_new<T> & B)
{
    if(A.nrows != X.nrows || A.ncols != X.ncols) eslog::error("matrix sizes dont match\n");
    if(B.nrows != X.nrows || B.ncols != X.ncols) eslog::error("matrix sizes dont match\n");
    if(A.order != X.order) eslog::error("matrix orders dont match\n");
    if(B.order != X.order) eslog::error("matrix orders dont match\n");

    T * vals_X = X.vals;
    T * vals_A = A.vals;
    T * vals_B = B.vals;
    size_t size_primary = X.get_size_primary();
    size_t size_secdary = X.get_size_secdary();
    bool move_start = ((X.uplo == 'U' && X.order == 'R') || (X.uplo == 'L' && X.order == 'C'));
    bool move_end   = ((X.uplo == 'L' && X.order == 'R') || (X.uplo == 'U' && X.order == 'C'));
    for(size_t i = 0; i < size_primary; i++) {
        size_t start = 0;
        size_t end = size_secdary;
        if(move_start) start = i;
        if(move_end) end = i;
        size_t size = end - start;
        T * sub_X = vals_X + i * X.ld;
        T * sub_A = vals_A + i * A.ld;
        T * sub_B = vals_B + i * B.ld;
        for(size_t j = start; j < end; j++) {
            sub_X[j] = alpha * sub_A[j] + beta * sub_B[j];
        }
    }
}



#define INSTANTIATE_T(T) \
template class lincomb_matrix_dnx<T>;

    #define INSTANTIATE \
    /* INSTANTIATE_T(float) */ \
    INSTANTIATE_T(double) \
    /* INSTANTIATE_T(std::complex<float>) */ \
    /* INSTANTIATE_T(std::complex<double>) */

        INSTANTIATE

    #undef INSTANTIATE
#undef INSTANTIATE_T



}
}
}
