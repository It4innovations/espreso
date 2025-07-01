
#ifndef SRC_MATH_OPERATIONS_MV_DNX_H
#define SRC_MATH_OPERATIONS_MV_DNX_H

#include "math/primitives_new/matrix_dense_view_new.h"
#include "math/primitives_new/vector_dense_view_new.h"



namespace espreso {
namespace math {
namespace operations {



template<typename T>
class mv_dnx
{
// y = alpha * A * x + beta * y
// respect symmetry and uplo of A
public:
    mv_dnx() = default;
    mv_dnx(const mv_dnx &) = delete;
    mv_dnx(mv_dnx &&) = default;
    mv_dnx & operator=(const mv_dnx &) = delete;
    mv_dnx & operator=(mv_dnx &&) = default;
    ~mv_dnx() = default;
public:
    void set_matrix_A(MatrixDenseView_new<T> * A_);
    void set_vector_x(VectorDenseView_new<T> * x_);
    void set_vector_y(VectorDenseView_new<T> * y_);
    void set_coefficients(T alpha_, T beta_);
    void perform();
    static void do_all(MatrixDenseView_new<T> * A, VectorDenseView_new<T> * x, VectorDenseView_new<T> * y, T alpha, T beta);
private:
    MatrixDenseView_new<T> * A = nullptr;
    VectorDenseView_new<T> * x = nullptr;
    VectorDenseView_new<T> * y = nullptr;
    T alpha = T{1};
    T beta = T{0};
};



}
}
}

#endif /* SRC_MATH_OPERATIONS_MV_DNX_H */
