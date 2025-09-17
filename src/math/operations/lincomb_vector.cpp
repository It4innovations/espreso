
#include "math/operations/lincomb_vector.h"

#include "math/operations/fill_dnx.h"
#include "basis/utilities/stacktimer.h"



namespace espreso {
namespace math {
namespace operations {



template<typename T>
void lincomb_vector<T>::set_vector_x(VectorDenseView_new<T> * x_)
{
    if(x != nullptr) eslog::error("vector x is already set\n");

    x = x_;
}



template<typename T>
void lincomb_vector<T>::set_vector_a(VectorDenseView_new<T> * a_)
{
    if(a != nullptr) eslog::error("vector a is already set\n");

    a = a_;
}



template<typename T>
void lincomb_vector<T>::set_vector_b(VectorDenseView_new<T> * b_)
{
    if(b != nullptr) eslog::error("vector b is already set\n");

    b = b_;
}



template<typename T>
void lincomb_vector<T>::set_coefficients(T alpha_, T beta_)
{
    alpha = alpha_;
    beta = beta_;
}



template<typename T>
void lincomb_vector<T>::perform()
{
    stacktimer::push("lincomb_vector::perform");

    if(x == nullptr) eslog::error("result vector x is not set\n");
    if(alpha != T{0} && a == nullptr) eslog::error("vector a is not set\n");
    if(beta != T{0} && b == nullptr) eslog::error("vector b is not set\n");
    if(!x->ator->is_data_accessible_cpu()) eslog::error("vector x must be cpu-accessible\n");
    if(a != nullptr && !a->ator->is_data_accessible_cpu()) eslog::error("vector a must be cpu-accessible\n");
    if(b != nullptr && !b->ator->is_data_accessible_cpu()) eslog::error("vector b must be cpu-accessible\n");

    if(alpha == T{0} && beta == T{0}) {
        lincomb_vector<T>::perform_zero(*x);
    }
    if(alpha != T{0} && beta == T{0}) {
        lincomb_vector<T>::perform_one(*x, alpha, *a);
    }
    if(alpha == T{0} && beta != T{0}) {
        lincomb_vector<T>::perform_one(*x, beta, *b);
    }
    if(alpha != T{0} && beta != T{0}) {
        lincomb_vector<T>::perform_two(*x, alpha, *a, beta, *b);
    }

    stacktimer::pop();
}



template<typename T>
void lincomb_vector<T>::do_all(VectorDenseView_new<T> * x, T alpha, VectorDenseView_new<T> * a, T beta, VectorDenseView_new<T> * b)
{
    lincomb_vector<T> instance;
    instance.set_vector_x(x);
    instance.set_vector_a(a);
    instance.set_vector_b(b);
    instance.set_coefficients(alpha, beta);
    instance.perform();
}



template<typename T>
void lincomb_vector<T>::perform_zero(VectorDenseView_new<T> & x)
{
    std::fill_n(x.vals, x.size, T{0});
}



template<typename T>
void lincomb_vector<T>::perform_one(VectorDenseView_new<T> & x, T alpha, VectorDenseView_new<T> & a)
{
    if(a.size != x.size) eslog::error("vector sizes dont match\n");

    T * vals_x = x.vals;
    T * vals_a = a.vals;

    for(size_t i = 0; i < x.size; i++) {
        vals_x[i] = alpha * vals_a[i];
    }
}



template<typename T>
void lincomb_vector<T>::perform_two(VectorDenseView_new<T> & x, T alpha, VectorDenseView_new<T> & a, T beta, VectorDenseView_new<T> & b)
{
    if(a.size != x.size) eslog::error("vector sizes dont match\n");
    if(b.size != x.size) eslog::error("vector sizes dont match\n");

    T * vals_x = x.vals;
    T * vals_a = a.vals;
    T * vals_b = b.vals;

    for(size_t i = 0; i < x.size; i++) {
        vals_x[i] = alpha * vals_a[i] + beta * vals_b[i];
    }
}



#define INSTANTIATE_T(T) \
template class lincomb_vector<T>;

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
