
#ifndef SRC_MATH_OPERATIONS_CONVERT_DNX_DNY_H
#define SRC_MATH_OPERATIONS_CONVERT_DNX_DNY_H

#include "math/primitives_new/matrix_dense_view_new.h"



namespace espreso {
namespace math {
namespace operations {



template<typename T>
class convert_dnx_dny
{
// does NOT respect uplo
public:
    convert_dnx_dny() = default;
    convert_dnx_dny(const convert_dnx_dny &) = delete;
    convert_dnx_dny(convert_dnx_dny &&) = default;
    convert_dnx_dny & operator=(const convert_dnx_dny &) = delete;
    convert_dnx_dny & operator=(convert_dnx_dny &&) = default;
    ~convert_dnx_dny() = default;
public:
    void set_matrix_src(MatrixDenseView_new<T> * M_src_);
    void set_matrix_dst(MatrixDenseView_new<T> * M_dst_);
    void set_conj(bool do_conj_);
    void perform();
    static void do_all(MatrixDenseView_new<T> * M_src, MatrixDenseView_new<T> * M_dst, bool do_conj);
private:
    MatrixDenseView_new<T> * M_src = nullptr;
    MatrixDenseView_new<T> * M_dst = nullptr;
    bool do_conj = true;
};



}
}
}

#endif /* SRC_MATH_OPERATIONS_CONVERT_DNX_DNY_H */
