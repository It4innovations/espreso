
#ifndef SRC_MATH_OPERATIONS_PIVOTS_TRAILS_CSX_H
#define SRC_MATH_OPERATIONS_PIVOTS_TRAILS_CSX_H

#include "math/primitives_new/matrix_csx_view_new.h"
#include "math/primitives_new/matrix_csx_data_new.h"
#include "math/primitives_new/vector_dense_view_new.h"



namespace espreso {
namespace math {
namespace operations {



template<typename T, typename I>
class pivots_trails_csx
{
public:
    void set_mode(char row_col_, char pivots_trails_, char completion_);
    void set_matrix(MatrixCsxView_new<T,I> * M_);
    void set_output_vector(VectorDenseView_new<I> * vec_);
    void perform();
    static void do_all(MatrixCsxView_new<T,I> * M, VectorDenseView_new<I> * vec, char row_col, char pivots_trails, char completion);
private:
    MatrixCsxView_new<T,I> * M = nullptr;
    VectorDenseView_new<I> * vec = nullptr;
    char row_col = '_';
    char pivots_trails = '_';
    char completion = 'N'; // None, Forward, Backward
private:
    static void perform_pivots(MatrixCsxView_new<T,I> & A, VectorDenseView_new<I> & pivots);
    static void perform_trails(MatrixCsxView_new<T,I> & A, VectorDenseView_new<I> & trails);
    static void complete_forward(VectorDenseView_new<I> & vec, I end_val);
    static void complete_backward(VectorDenseView_new<I> & vec, I end_val);
};



}
}
}

#endif /* SRC_MATH_OPERATIONS_PIVOTS_TRAILS_CSX_H */
