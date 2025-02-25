
#ifndef SRC_MATH_OPERATIONS_AUXILIARY_SORT_CSX_CSX_PIVTRL
#define SRC_MATH_OPERATIONS_AUXILIARY_SORT_CSX_CSX_PIVTRL

#include "math/primitives_new/matrix_csx_view_new.h"



namespace espreso {
namespace math {
namespace operations {



template<typename T, typename I>
class sort_csx_csx_pivtrl
{
public:
    void set_matrix_src(MatrixCsxView_new<T,I> * M_src_);
    void set_matrix_dst(MatrixCsxView_new<T,I> * M_dst_);
    void set_mode(char row_col_, char piv_trl_, char asc_dsc_);
    void perform();
    static void do_all(MatrixCsxView_new<T,I> * M_src, MatrixCsxView_new<T,I> * M_dst, char row_col, char piv_trl, char asc_dsc);
private:
    MatrixCsxView_new<T,I> * M_src = nullptr;
    MatrixCsxView_new<T,I> * M_dst = nullptr;
    char row_col = '_';
    char piv_trl = '_';
    char asc_dsc = '_';
    bool mode_set = false;
};



}
}
}

#endif /* SRC_MATH_OPERATIONS_AUXILIARY_SORT_CSX_CSX_PIVTRL */
