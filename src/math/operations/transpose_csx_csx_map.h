
#ifndef SRC_MATH_OPERATIONS_TRANSPOSE_CSX_CSX_MAP_H
#define SRC_MATH_OPERATIONS_TRANSPOSE_CSX_CSX_MAP_H

#include "math/primitives_new.h"



namespace espreso {
namespace math {
namespace operations {



template<typename T, typename I>
class transpose_csx_csx_map
{
    // ignores uplo
public:
    transpose_csx_csx_map() = default;
    transpose_csx_csx_map(const transpose_csx_csx_map &) = delete;
    transpose_csx_csx_map(transpose_csx_csx_map &&) = default;
    transpose_csx_csx_map & operator=(const transpose_csx_csx_map &) = delete;
    transpose_csx_csx_map & operator=(transpose_csx_csx_map &&) = default;
    ~transpose_csx_csx_map() = default;
public:
    void set_matrix_src(MatrixCsxView_new<T,I> * M_src_);
    void set_matrix_dst(MatrixCsxView_new<T,I> * M_dst_);
    void perform_pattern();
    void perform_values();
    void perform_all();
    static void do_all(MatrixCsxView_new<T,I> * M_src, MatrixCsxView_new<T,I> * M_dst);
private:
    MatrixCsxView_new<T,I> * M_src = nullptr;
    MatrixCsxView_new<T,I> * M_dst = nullptr;
    VectorDenseData_new<I> map;
    bool perform_pattern_called = false;
};



}
}
}

#endif /* SRC_MATH_OPERATIONS_TRANSPOSE_CSX_CSX_MAP_H */
