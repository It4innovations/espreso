
#ifndef SRC_MATH_OPERATIONS_LINCOMB_VECTOR_H
#define SRC_MATH_OPERATIONS_LINCOMB_VECTOR_H

#include "math/primitives_new/vector_dense_view_new.h"



namespace espreso {
namespace math {
namespace operations {



template<typename T>
class lincomb_vector
{
// x = alpha * a + beta * b
// inplace allowed (x=a, x=b)
public:
    lincomb_vector() = default;
    lincomb_vector(const lincomb_vector &) = delete;
    lincomb_vector(lincomb_vector &&) = default;
    lincomb_vector & operator=(const lincomb_vector &) = delete;
    lincomb_vector & operator=(lincomb_vector &&) = default;
    ~lincomb_vector() = default;
public:
    void set_vector_x(VectorDenseView_new<T> * x_);
    void set_vector_a(VectorDenseView_new<T> * a_);
    void set_vector_b(VectorDenseView_new<T> * b_);
    void set_coefficients(T alpha_, T beta_);
    void perform();
    static void do_all(VectorDenseView_new<T> * x, T alpha, VectorDenseView_new<T> * a, T beta, VectorDenseView_new<T> * b);
private:
    VectorDenseView_new<T> * x = nullptr;
    VectorDenseView_new<T> * a = nullptr;
    VectorDenseView_new<T> * b = nullptr;
    T alpha = T{0};
    T beta = T{0};
private:
    static void perform_zero(VectorDenseView_new<T> & x);
    static void perform_one(VectorDenseView_new<T> & x, T alpha, VectorDenseView_new<T> & a);
    static void perform_two(VectorDenseView_new<T> & x, T alpha, VectorDenseView_new<T> & a, T beta, VectorDenseView_new<T> & b);
};



}
}
}

#endif /* SRC_MATH_OPERATIONS_LINCOMB_VECTOR_H */
