
#include "math/operations/quadrisect_csx_csy.h"

#include "basis/utilities/stacktimer.h"



namespace espreso {
namespace math {
namespace operations {



template<typename T, typename I>
void quadrisect_csx_csy<T,I>::set_matrix_src(MatrixCsxView_new<T,I> * M_src_)
{
    if(M_src != nullptr) eslog::error("matrix M_src is already set\n");

    M_src = M_src_;
}



template<typename T, typename I>
void quadrisect_csx_csy<T,I>::set_bounds(size_t row_cut_, size_t col_cut_)
{
    if(called_set_bounds) eslog::error("bounds are already set\n");

    row_cut = row_cut_;
    col_cut = col_cut_;

    called_set_bounds = true;
}



template<typename T, typename I>
void quadrisect_csx_csy<T,I>::setup()
{
    stacktimer::push("quadrisect_csx_csy::setup");

    if(called_setup) eslog::error("setup has already been called\n");
    if(!called_set_bounds) eslog::error("bounds are not set\n");
    if(M_src == nullptr) eslog::error("source matrix is not set\n");

    if(row_cut > M_src->nrows || col_cut > M_src->ncols) eslog::error("wrong bounds\n");

    op_sub_11.set_matrix_src(M_src);
    op_sub_11.set_bounds(0, row_cut, 0, col_cut);
    op_sub_11.setup();

    op_sub_12.set_matrix_src(M_src);
    op_sub_12.set_bounds(0, row_cut, col_cut, M_src->ncols);
    op_sub_12.setup();

    op_sub_21.set_matrix_src(M_src);
    op_sub_21.set_bounds(row_cut, M_src->nrows, 0, col_cut);
    op_sub_21.setup();

    op_sub_22.set_matrix_src(M_src);
    op_sub_22.set_bounds(row_cut, M_src->nrows, col_cut, M_src->ncols);
    op_sub_22.setup();

    stacktimer::pop();

    called_setup = true;
}



template<typename T, typename I>
size_t quadrisect_csx_csy<T,I>::get_output_matrix_11_nnz()
{
    if(!called_setup) eslog::error("setup was not called\n");

    return op_sub_11.get_output_matrix_nnz();
}



template<typename T, typename I>
size_t quadrisect_csx_csy<T,I>::get_output_matrix_12_nnz()
{
    if(!called_setup) eslog::error("setup was not called\n");

    return op_sub_12.get_output_matrix_nnz();
}



template<typename T, typename I>
size_t quadrisect_csx_csy<T,I>::get_output_matrix_21_nnz()
{
    if(!called_setup) eslog::error("setup was not called\n");

    return op_sub_21.get_output_matrix_nnz();
}



template<typename T, typename I>
size_t quadrisect_csx_csy<T,I>::get_output_matrix_22_nnz()
{
    if(!called_setup) eslog::error("setup was not called\n");

    return op_sub_22.get_output_matrix_nnz();
}



template<typename T, typename I>
void quadrisect_csx_csy<T,I>::set_matrices_dst(MatrixCsxView_new<T,I> * M_dst_11_, MatrixCsxView_new<T,I> * M_dst_12_, MatrixCsxView_new<T,I> * M_dst_21_, MatrixCsxView_new<T,I> * M_dst_22_)
{
    if(called_set_dst) eslog::error("dst matrices are already set\n");

    M_dst_11 = M_dst_11_;
    M_dst_12 = M_dst_12_;
    M_dst_21 = M_dst_21_;
    M_dst_22 = M_dst_22_;

    if(M_dst_11 != nullptr && (M_dst_11->nrows != row_cut || M_dst_11->ncols != col_cut || M_dst_11->nnz != op_sub_11.get_output_matrix_nnz())) eslog::error("wrong dst 11 matrix size\n");
    if(M_dst_12 != nullptr && (M_dst_12->nrows != row_cut || M_dst_12->ncols != M_src->ncols - col_cut || M_dst_12->nnz != op_sub_12.get_output_matrix_nnz())) eslog::error("wrong dst 12 matrix size\n");
    if(M_dst_21 != nullptr && (M_dst_21->nrows != M_src->nrows - row_cut || M_dst_21->ncols != col_cut || M_dst_21->nnz != op_sub_21.get_output_matrix_nnz())) eslog::error("wrong dst 21 matrix size\n");
    if(M_dst_22 != nullptr && (M_dst_22->nrows != M_src->nrows - row_cut || M_dst_22->ncols != M_src->ncols - col_cut || M_dst_22->nnz != op_sub_22.get_output_matrix_nnz())) eslog::error("wrong dst 22 matrix size\n");

    op_sub_11.set_matrix_dst(M_dst_11);
    op_sub_12.set_matrix_dst(M_dst_12);
    op_sub_21.set_matrix_dst(M_dst_21);
    op_sub_22.set_matrix_dst(M_dst_22);

    called_set_dst = true;
}



template<typename T, typename I>
void quadrisect_csx_csy<T,I>::perform()
{
    stacktimer::push("quadrisect_csx_csy::perform");

    if(!called_setup) eslog::error("setup was not called\n");
    if(!called_set_dst) eslog::error("destination matrices are not set\n");

    if(M_dst_11 != nullptr) op_sub_11.perform();
    if(M_dst_12 != nullptr) op_sub_12.perform();
    if(M_dst_21 != nullptr) op_sub_21.perform();
    if(M_dst_22 != nullptr) op_sub_22.perform();

    stacktimer::pop();
}



template<typename T, typename I>
void quadrisect_csx_csy<T,I>::do_all(MatrixCsxView_new<T,I> * M_src, MatrixCsxView_new<T,I> * M_dst_11, MatrixCsxView_new<T,I> * M_dst_12, MatrixCsxView_new<T,I> * M_dst_21, MatrixCsxView_new<T,I> * M_dst_22, size_t row_cut, size_t col_cut)
{
    quadrisect_csx_csy<T,I> instance;
    instance.set_matrix_src(M_src);
    instance.set_matrices_dst(M_dst_11, M_dst_12, M_dst_21, M_dst_22);
    instance.set_bounds(row_cut, col_cut);
    instance.setup();
    instance.perform();
}



#define INSTANTIATE_T_I(T,I) \
template class quadrisect_csx_csy<T,I>;

    #define INSTANTIATE_T(T) \
    INSTANTIATE_T_I(T,int32_t) \
    /* INSTANTIATE_T_I(T,int64_t) */

        #define INSTANTIATE \
        /* INSTANTIATE_T(float) */ \
        INSTANTIATE_T(double) \
        /* INSTANTIATE_T(std::complex<float>) */ \
        INSTANTIATE_T(std::complex<double>)

            INSTANTIATE

        #undef INSTANTIATE
    #undef INSTANTIATE_T
#undef INSTANTIATE_T_I



}
}
}
