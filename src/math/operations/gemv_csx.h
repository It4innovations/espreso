
#ifndef SRC_MATH_OPERATIONS_GEMV_DNX_H
#define SRC_MATH_OPERATIONS_GEMV_DNX_H

#include "math/primitives_new.h"



namespace espreso {
namespace math {
namespace operations {



template<typename T, typename I>
class gemv_csx
{
    // ignores uplo and symmetry
    // y = alpha * A * x + beta * y
public:
    gemv_csx() = default;
    gemv_csx(const gemv_csx &) = delete;
    gemv_csx(gemv_csx &&) = default;
    gemv_csx & operator=(const gemv_csx &) = delete;
    gemv_csx & operator=(gemv_csx &&) = default;
    ~gemv_csx() = default;
public:
    void set_matrix_A(MatrixCsxView_new<T,I> * A_);
    void set_vector_x(VectorDenseView_new<T> * x_);
    void set_vector_y(VectorDenseView_new<T> * y_);
    void set_coefficients(T alpha_, T beta_);
    void perform();
    static void do_all(MatrixCsxView_new<T,I> * A, VectorDenseView_new<T> * x, VectorDenseView_new<T> * y, T alpha, T beta);
private:
    MatrixCsxView_new<T,I> * A = nullptr;
    VectorDenseView_new<T> * x = nullptr;
    VectorDenseView_new<T> * y = nullptr;
    T alpha = T{1};
    T beta = T{0};
};



}
}
}



#endif /* SRC_MATH_OPERATIONS_GEMV_DNX_H */
