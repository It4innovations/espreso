

#include "math/operations/gemv_csx.h"

#include "basis/utilities/stacktimer.h"



namespace espreso {
namespace math {
namespace operations {



template<typename T, typename I>
void gemv_csx<T,I>::set_matrix_A(MatrixCsxView_new<T,I> * A_)
{
    if(A != nullptr) eslog::error("matrix A is already set\n");

    A = A_;
}



template<typename T, typename I>
void gemv_csx<T,I>::set_vector_x(VectorDenseView_new<T> * x_)
{
    if(x != nullptr) eslog::error("vector x is already set\n");

    x = x_;
}



template<typename T, typename I>
void gemv_csx<T,I>::set_vector_y(VectorDenseView_new<T> * y_)
{
    if(y != nullptr) eslog::error("vector y is already set\n");

    y = y_;
}



template<typename T, typename I>
void gemv_csx<T,I>::set_coefficients(T alpha_, T beta_)
{
    alpha = alpha_;
    beta = beta_;
}



template<typename T, typename I>
void gemv_csx<T,I>::perform()
{
    stacktimer::push("gemv_csx::perform");

    if(A == nullptr) eslog::error("matrix A is not set\n");
    if(x == nullptr) eslog::error("vector x is not set\n");
    if(y == nullptr) eslog::error("vector y is not set\n");
    if(!A->ator->is_data_accessible_cpu()) eslog::error("matrix A must be cpu-accessible\n");
    if(!x->ator->is_data_accessible_cpu()) eslog::error("vector x must be cpu-accessible\n");
    if(!y->ator->is_data_accessible_cpu()) eslog::error("vector y must be cpu-accessible\n");
    if(x->size != A->ncols) eslog::error("wrong x size\n");
    if(y->size != A->nrows) eslog::error("wrong y size\n");

    // TODO: use sparse blas libraries instead

    T * x_vals = x->vals;
    T * y_vals = y->vals;

    if(A->order == 'R') {
        if(beta == T{0}) {
            std::fill_n(y_vals, y->size, T{0});
        }
        I nrows = A->nrows;
        for(I row = 0; row < nrows; row++) {
            I start = A->ptrs[row];
            I end = A->ptrs[row+1];
            T y_val = T{0};
            for(I i = start; i < end; i++) {
                I col = A->idxs[i];
                T val = A->vals[i];
                y_val += val * x_vals[col];
            }
            y_vals[row] = alpha * y_val + beta * y_vals[row];
        }
    }
    else {
        if(beta == T{0}) {
            std::fill_n(y_vals, y->size, T{0});
        }
        else {
            for(size_t i = 0; i < y->size; i++) {
                y_vals[i] *= beta;
            }
        }
        I ncols = A->ncols;
        for(I col = 0; col < ncols; col++) {
            I start = A->ptrs[col];
            I end = A->ptrs[col+1];
            T x_val = x_vals[col];
            for(I i = start; i < end; i++) {
                I row = A->idxs[i];
                T val = A->vals[i];
                y_vals[row] += alpha * val * x_val;
            }
        }
    }

    stacktimer::pop();
}



template<typename T, typename I>
void gemv_csx<T,I>::do_all(MatrixCsxView_new<T,I> * A, VectorDenseView_new<T> * x, VectorDenseView_new<T> * y, T alpha, T beta)
{
    gemv_csx<T,I> instance;
    instance.set_matrix_A(A);
    instance.set_vector_x(x);
    instance.set_vector_y(y);
    instance.set_coefficients(alpha, beta);
    instance.perform();
}



#define INSTANTIATE_T_I(T,I) \
template class gemv_csx<T,I>;

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
