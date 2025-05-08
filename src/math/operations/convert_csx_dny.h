
#ifndef SRC_MATH_OPERATIONS_CONVERT_CSX_DNY_H
#define SRC_MATH_OPERATIONS_CONVERT_CSX_DNY_H

#include "math/primitives_new/matrix_csx_view_new.h"
#include "math/primitives_new/matrix_dense_view_new.h"



namespace espreso {
namespace math {
namespace operations {



template<typename T, typename I>
class convert_csx_dny
{
// respect uplo
public:
    convert_csx_dny() = default;
    convert_csx_dny(const convert_csx_dny &) = delete;
    convert_csx_dny(convert_csx_dny &&) = default;
    convert_csx_dny & operator=(const convert_csx_dny &) = delete;
    convert_csx_dny & operator=(convert_csx_dny &&) = default;
    ~convert_csx_dny() = default;
public:
    void set_matrix_src(MatrixCsxView_new<T,I> * M_src_);
    void set_matrix_dst(MatrixDenseView_new<T> * M_dst_);
    void perform_zerofill();
    void perform_copyvals();
    void perform_all();
    static void do_all(MatrixCsxView_new<T,I> * M_src, MatrixDenseView_new<T> * M_dst);
private:
    MatrixCsxView_new<T,I> * M_src = nullptr;
    MatrixDenseView_new<T> * M_dst = nullptr;
    bool zerofill_called = false;
};



}
}
}

#endif /* SRC_MATH_OPERATIONS_CONVERT_CSX_DNY_H */
