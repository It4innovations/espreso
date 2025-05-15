
#include "math/operations/submatrix_dnx_dnx_view.h"

#include "math/primitives_new/matrix_dense_data_new.h"
#include "basis/utilities/stacktimer.h"



namespace espreso {
namespace math {
namespace operations {



template<typename T>
void submatrix_dnx_dnx_view<T>::set_matrix_src(MatrixDenseView_new<T> * M_src_)
{
    if(M_src != nullptr) eslog::error("matrix M_src is already set\n");

    M_src = M_src_;
}



template<typename T>
void submatrix_dnx_dnx_view<T>::set_matrix_dst(MatrixDenseView_new<T> * M_dst_)
{
    if(M_dst != nullptr) eslog::error("matrix M_dst is already set\n");

    M_dst = M_dst_;
}



template<typename T>
void submatrix_dnx_dnx_view<T>::set_bounds(size_t row_start_, size_t row_end_, size_t col_start_, size_t col_end_)
{
    row_start = row_start_;
    row_end = row_end_;
    col_start = col_start_;
    col_end = col_end_;
    num_rows = row_end - row_start;
    num_cols = col_end - col_start;

    called_set_bounds = true;
}



template<typename T>
void submatrix_dnx_dnx_view<T>::perform()
{
    stacktimer::push("submatrix_dnx_dnx_view::perform");

    if(M_src == nullptr) eslog::error("source matrix has not been set\n");
    if(M_dst == nullptr) eslog::error("destination matrix has not been set\n");
    if(M_dst->nrows != num_rows || M_dst->ncols != num_cols) eslog::error("destination matrix size does not match bounds\n");
    if(!called_set_bounds) eslog::error("bounds are not set\n");
    if(row_start > row_end || row_end > M_src->nrows || col_start > col_end || col_end > M_src->ncols) eslog::error("wrong bounds\n");
    if(M_src->order != M_dst->order) eslog::error("ordes must match\n");
    if(M_src->ld != M_dst->ld) eslog::error("leading dimensions must match\n");
    if(dynamic_cast<MatrixDenseData_new<T>*>(M_dst) != nullptr) eslog::error("cannot be used with data as destination\n");

    M_dst->vals = M_src->vals + row_start * M_src->get_stride_row() + col_start * M_dst->get_stride_col();

    stacktimer::pop();
}



template<typename T>
void submatrix_dnx_dnx_view<T>::do_all(MatrixDenseView_new<T> * M_src, MatrixDenseView_new<T> * M_dst, size_t row_start, size_t row_end, size_t col_start, size_t col_end)
{
    submatrix_dnx_dnx_view<T> instance;
    instance.set_matrix_src(M_src);
    instance.set_matrix_dst(M_dst);
    instance.set_bounds(row_start, row_end, col_start, col_end);
    instance.perform();
}



#define INSTANTIATE_T(T) \
template class submatrix_dnx_dnx_view<T>;

    #define INSTANTIATE \
    /* INSTANTIATE_T(float) */ \
    INSTANTIATE_T(double) \
    /* INSTANTIATE_T(std::complex<float>) */ \
    INSTANTIATE_T(std::complex<double>)

        INSTANTIATE

    #undef INSTANTIATE
#undef INSTANTIATE_T



}
}
}
