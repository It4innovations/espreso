
#ifndef SRC_MATH_OPERATIONS_COVNERT_CSX_CSY_MAP_H
#define SRC_MATH_OPERATIONS_COVNERT_CSX_CSY_MAP_H

#include "math/primitives_new/matrix_csx_view_new.h"
#include "math/primitives_new/vector_dense_data_new.h"



namespace espreso {
namespace math {
namespace operations {



template<typename T, typename I>
class convert_csx_csy_map
{
public:
    convert_csx_csy_map() = default;
    convert_csx_csy_map(const convert_csx_csy_map &) = delete;
    convert_csx_csy_map(convert_csx_csy_map &&) = delete;
    convert_csx_csy_map & operator=(const convert_csx_csy_map &) = delete;
    convert_csx_csy_map & operator=(convert_csx_csy_map &&) = delete;
    ~convert_csx_csy_map();
public:
    void set_matrix_src(MatrixCsxView_new<T,I> * M_src_);
    void set_matrix_dst(MatrixCsxView_new<T,I> * M_dst_);
    void perform_pattern();
    void perform_values();
    void perform_all();
    void finalize();
    static void do_all(MatrixCsxView_new<T,I> * M_src, MatrixCsxView_new<T,I> * M_dst);
private:
    MatrixCsxView_new<T,I> * M_src = nullptr;
    MatrixCsxView_new<T,I> * M_dst = nullptr;
    VectorDenseData<I> map;
    bool perform_pattern_called = false;
};



}
}
}

#endif /* SRC_MATH_OPERATIONS_COVNERT_CSX_CSY_MAP_H */
