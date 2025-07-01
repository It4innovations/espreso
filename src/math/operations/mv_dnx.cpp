
#include "math/operations/mv_dnx.h"

#include "math/wrappers/math.blas.h"
#include "basis/utilities/stacktimer.h"



namespace espreso {
namespace math {
namespace operations {



template<typename T>
void mv_dnx<T>::set_matrix_A(MatrixDenseView_new<T> * A_)
{
    if(A != nullptr) eslog::error("matrix A is already set\n");

    A = A_;
}



template<typename T>
void mv_dnx<T>::set_vector_x(VectorDenseView_new<T> * x_)
{
    if(x != nullptr) eslog::error("vector x is already set\n");

    x = x_;
}



template<typename T>
void mv_dnx<T>::set_vector_y(VectorDenseView_new<T> * y_)
{
    if(y != nullptr) eslog::error("vector y is already set\n");

    y = y_;
}



template<typename T>
void mv_dnx<T>::set_coefficients(T alpha_, T beta_)
{
    alpha = alpha_;
    beta = beta_;
}



template<typename T>
void mv_dnx<T>::perform()
{
    stacktimer::push("mv_dnx::perform");

    if(A == nullptr) eslog::error("matrix A is not set\n");
    if(x == nullptr) eslog::error("vector x is not set\n");
    if(y == nullptr) eslog::error("vector y is not set\n");
    if(x->size != A->ncols) eslog::error("incompatible vector x\n");
    if(y->size != A->nrows) eslog::error("incompatible vector y\n");

    if(is_hermitian<T>(A->prop.symm)) {
        blas::apply_hermitian<T,int>(*y, alpha, *A, beta, *x);
    }
    else if(is_symmetric<T>(A->prop.symm)) {
        eslog::error("symmetric matrices not supported yet\n");
        // but real symmetric falls into hermitian so in most cases this is not a problem
    }
    else {
        blas::apply<T,int>(*y, alpha, *A, beta, *x);
    }

    stacktimer::pop();
}



template<typename T>
void mv_dnx<T>::do_all(MatrixDenseView_new<T> * A, VectorDenseView_new<T> * x, VectorDenseView_new<T> * y, T alpha, T beta)
{
    mv_dnx<T> instance;
    instance.set_matrix_A(A);
    instance.set_vector_x(x);
    instance.set_vector_y(y);
    instance.set_coefficients(alpha, beta);
    instance.perform();
}



#define INSTANTIATE_T(T) \
template class mv_dnx<T>;

    #define INSTANTIATE \
    /* INSTANTIATE_T(float) */ \
    INSTANTIATE_T(double) \
    /* INSTANTIATE_T(std::complex<float>) */ \
    INSTANTIATE_T(std::complex<double>)

        INSTANTIATE

    #undef INSTANTIATE
#undef INSTANTIATE_T



}
}
}
