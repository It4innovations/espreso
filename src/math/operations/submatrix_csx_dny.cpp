
#include "math/operations/submatrix_csx_dny.h"



namespace espreso {
namespace math {
namespace operations {



template<typename T, typename I>
void submatrix_csx_dny<T,I>::set_matrix_src(const MatrixCsxView_new<T,I> & M_src_)
{
    M_src = M_src_;
}



template<typename T, typename I>
void submatrix_csx_dny<T,I>::set_matrix_dst(const MatrixDenseView_new<T> & M_dst_)
{
    M_dst = M_dst_;
}



template<typename T, typename I>
void submatrix_csx_dny<T,I>::set_bounds(size_t row_start_, size_t row_end_, size_t col_start_, size_t col_end_)
{
    if(M_src == nullptr) eslog::error("source matrix is not set\n");

    row_start = row_start_;
    row_end = row_end_;
    col_start = col_start_;
    col_end = col_end_;
    num_rows = end_row - start_row;
    num_cols = end_col - start_col;

    if(row_start > row_end || row_end > M_src->nrows || col_start > col_end || col_end > M_src.ncols) eslog::error("wrong bounds\n");

    bound_set = true;
}



template<typename T, typename I>
void submatrix_csx_dny<T,I>::perform_zerofill()
{
    if(M_src == nullptr) eslog::error("source matrix is not set\n");
    if(M_dst == nullptr) eslog::error("destination matrix is not set\n");
    if(!bound_set) eslog::error("bounds are not set\n");
    if(M_dst->nrows != num_rows || M_dst->ncols != num_cols) eslog::error("wrong output matrix size\n");
    
    size_t num_blocks = M_dst->get_num_blocks();
    size_t block_size = M_dst->get_block_size();
    for(size_t i = 0; i < num_blocks; i++) {
        std::fill_n(M_dst->vals + i * M_dst->ld, block_size, T{0});
    }

    zerofill_called = true;
}



template<typename T, typename I>
void submatrix_csx_dny<T,I>::perform_copyvals()
{
    if(M_src == nullptr) eslog::error("source matrix is not set\n");
    if(M_dst == nullptr) eslog::error("destination matrix is not set\n");
    if(!bound_set) eslog::error("bounds are not set\n");
    if(!zerofill_called) eslog::error("zerofill was not called\n");
    if(M_dst->nrows != num_rows || M_dst->ncols != num_cols) eslog::error("wrong output matrix size\n");

    size_t start_prim = 0;
    size_t end_prim = 0;
    size_t start_sec = 0;
    size_t end_sec = 0;
    if(M_src->order == 'R') {
        start_prim = start_row;
        end_prim = end_row;
        start_sec = start_col;
        end_sec = end_col;
    }
    if(M_src->order == 'C') {
        start_prim = start_col;
        end_prim = end_col;
        start_sec = start_row;
        end_sec = end_row;
    }

    size_t stride_prim = ((M_src->order == M_dst->order) ? M_dst->ld : 1);
    size_t stride_sec = ((M_src->order == M_dst->order) ? 1 : M_dst->ld);

    I * srcptrs = M_src->ptrs;
    I * srcidxs = M_src->idxs;
    T * srcvals = M_src->vals;
    T * dstvals = M_dst->vals;

    for(size_t ip = start_prim; ip < end_prim; ip++)
    {
        I start = srcptrs[ip];
        I end = srcptrs[ip+1];
        I i = start;
        while(i < end && srcidxs[i] < start_sec) {
            i++;
        }
        while(i < end && srcidxs[i] < end_sec) {
            I is = srcidxs[i];
            T v = srcvals[i];
            dstvals[ip * stride_prim + is * stride_sec] = v;
            i++;
        }
    }
}



template<typename T, typename I>
void submatrix_csx_dny<T,I>::perform_all()
{
    perform_zerofill();
    perform_copyvals();
}



template<typename T, typename I>
void submatrix_csx_dny<T,I>::do_all(MatrixCsxView_new<T,I> * M_src, MatrixDenseView_new<T> * M_dst, size_t row_start, size_t row_end, size_t col_start, size_t col_end)
{
    submatrix_csx_dny<T,I> instance;
    instance.set_matrix_src(M_src);
    instance.set_matrix_dst(M_dst);
    instance.set_bounds(row_start, row_end, col_start, col_end);
    instance.perform_zerofill();
    instance.perform_copyvals();
}



#define INSTANTIATE_T_I(T,I) \
template class submatrix_csx_dny<T,I>;

    #define INSTANTIATE_T(T) \
    INSTANTIATE_T_I(T,int32_t) \
    /* INSTANTIATE_T_I(T,int64_t) */

        #define INSTANTIATE \
        /* INSTANTIATE_T(float) */ \
        INSTANTIATE_T(double) \
        /* INSTANTIATE_T(std::complex<float>) */ \
        /* INSTANTIATE_T(std::complex<double>) */

            INSTANTIATE

        #undef INSTANTIATE
    #undef INSTANTIATE_T
#undef INSTANTIATE_T_I



}
}
}
