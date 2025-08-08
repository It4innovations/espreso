
#ifndef SRC_MATH_OPERATIONS_COMPLETE_CSX_CSX_H
#define SRC_MATH_OPERATIONS_COMPLETE_CSX_CSX_H

#include "math/primitives_new.h"



namespace espreso {
namespace math {
namespace operations {



template<typename T, typename I>
class complete_csx_csx_map
{
public:
    complete_csx_csx_map() = default;
    complete_csx_csx_map(const complete_csx_csx_map &) = delete;
    complete_csx_csx_map(complete_csx_csx_map &&) = default;
    complete_csx_csx_map & operator=(const complete_csx_csx_map &) = delete;
    complete_csx_csx_map & operator=(complete_csx_csx_map &&) = default;
    ~complete_csx_csx_map() = default;
public:
    void set_matrix_src(MatrixCsxView_new<T,I> * M_src_);
    void set_matrix_dst(MatrixCsxView_new<T,I> * M_dst_);
    void set_conj(bool do_conj_);
    size_t get_dst_nnz();
    void perform_pattern();
    void perform_values();
private:
    MatrixCsxView_new<T,I> * M_src = nullptr;
    MatrixCsxView_new<T,I> * M_dst = nullptr;
    VectorDenseData_new<I> map;
    bool do_conj = true;
    bool called_perform_pattern = false;
};



}
}
}



#endif /* SRC_MATH_OPERATIONS_COMPLETE_CSX_CSX_H */
