
#ifndef SRC_MATH_OPERATIONS_AUXILIARY_PIVOTS_TRAILS_CSX_H
#define SRC_MATH_OPERATIONS_AUXILIARY_PIVOTS_TRAILS_CSX_H

#include "math/primitives_new/matrix_csx_view_new.h"
#include "math/primitives_new/matrix_csx_data_new.h"
#include "math/primitives_new/vector_dense_view_new.h"



template<typename T, typename I>
class pivots_trails_csx
{
public:
    void set_mode(char row_col_, char pivots_trails_);
    void set_matrix(MatrixCsxView_new<T,I> * M_);
    void set_output_vector(VectorDenseView_new<I> * vec_);
    void perform();
    static void do_all(MatrixCsxView_new<T,I> * M, VectorDenseView_new<I> * vec, char row_col, char pivots_trails);
private:
    MatrixCsxView_new<T,I> * M = nullptr;
    VectorDenseView_new<I> * vec = nullptr;
    char row_col = '_';
    char pivots_trails = '_';
private:
    void perform_pivots(MatrixCsxView_new<T,I> & A, VectorDenseView_new<I> & pivots);
    void perform_trails(MatrixCsxView_new<T,I> & A, VectorDenseView_new<I> & trails);
};



#endif /* SRC_MATH_OPERATIONS_AUXILIARY_PIVOTS_TRAILS_CSX_H */
