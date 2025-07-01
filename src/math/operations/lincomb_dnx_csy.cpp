
#include "math/operations/lincomb_dnx_csy.h"

#include "math/operations/lincomb_matrix_dnx.h"
#include "basis/utilities/stacktimer.h"



namespace espreso {
namespace math {
namespace operations {



template<typename T, typename I>
void lincomb_dnx_csy<T,I>::set_matrix_X(MatrixDenseView_new<T> * X_)
{
    if(X != nullptr) eslog::error("matrix X is already set\n");

    X = X_;
}



template<typename T, typename I>
void lincomb_dnx_csy<T,I>::set_matrix_A(MatrixDenseView_new<T> * A_)
{
    if(A != nullptr) eslog::error("matrix A is already set\n");

    A = A_;
}



template<typename T, typename I>
void lincomb_dnx_csy<T,I>::set_matrix_B(MatrixCsxView_new<T,I> * B_)
{
    if(B != nullptr) eslog::error("matrix B is already set\n");

    B = B_;
}



template<typename T, typename I>
void lincomb_dnx_csy<T,I>::set_coefficients(T alpha_, T beta_)
{
    alpha = alpha_;
    beta = beta_;
}



template<typename T, typename I>
void lincomb_dnx_csy<T,I>::perform()
{
    stacktimer::push("lincomb_dnx_csy::perform");

    if(X == nullptr) eslog::error("result matrix X is not set\n");
    if(alpha != T{0} && A == nullptr) eslog::error("matrix A is not set\n");
    if(beta != T{0} && B == nullptr) eslog::error("matrix B is not set\n");
    if(!X->ator->is_data_accessible_cpu()) eslog::error("matrix X must be cpu-accessible\n");
    if(A != nullptr && !A->ator->is_data_accessible_cpu()) eslog::error("matrix A must be cpu-accessible\n");
    if(B != nullptr && !B->ator->is_data_accessible_cpu()) eslog::error("matrix B must be cpu-accessible\n");
    if(A != nullptr && (A->nrows != X->nrows || A->ncols != X->ncols)) eslog::error("nonmatching size of matrix A\n");
    if(B != nullptr && (B->nrows != X->nrows || B->ncols != X->ncols)) eslog::error("nonmatching size of matrix B\n");
    if(A != nullptr && A->order != X->order) eslog::error("order of matrices A and X must match\n");
    if(A != nullptr && A->prop.uplo != X->prop.uplo) eslog::error("A matrix uplo does not match X\n");
    if(B != nullptr && B->prop.uplo != X->prop.uplo) eslog::error("B matrix uplo does not match X\n");

    if(!(X == A && alpha == T{1})) {
        lincomb_matrix_dnx<T>::do_all(X, alpha, A, 0, nullptr);
    }

    if(beta != T{0}) {
        size_t primary_size = B->get_size_primary();
        size_t dstld = X->ld;
        I * srcptrs = B->ptrs;
        I * srcidxs = B->idxs;
        T * srcvals = B->vals;
        T * dstvals = X->vals;

        if(B->order == X->order)
        {
            for(size_t r = 0; r < primary_size; r++)
            {
                I start = srcptrs[r];
                I end = srcptrs[r+1];
                T * row = dstvals + r * dstld;
                for(I i = start; i < end; i++)
                {
                    I c = srcidxs[i];
                    T v = srcvals[i];
                    row[c] += beta * v;
                }
            }
        }
        else
        {
            for(size_t r = 0; r < primary_size; r++)
            {
                I start = srcptrs[r];
                I end = srcptrs[r+1];
                for(I i = start; i < end; i++)
                {
                    I c = srcidxs[i];
                    T v = srcvals[i];
                    dstvals[r + dstld * c] += beta * v;
                }
            }
        }
    }

    stacktimer::pop();
}



template<typename T, typename I>
void lincomb_dnx_csy<T,I>::do_all(MatrixDenseView_new<T> * X, T alpha, MatrixDenseView_new<T> * A, T beta, MatrixCsxView_new<T,I> * B)
{
    lincomb_dnx_csy<T,I> instance;
    instance.set_matrix_X(X);
    instance.set_matrix_A(A);
    instance.set_matrix_B(B);
    instance.set_coefficients(alpha, beta);
    instance.perform();
}



#define INSTANTIATE_T_I(T,I) \
template class lincomb_dnx_csy<T,I>;

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
