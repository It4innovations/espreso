
#include "math/operations/submatrix_csx_csx.h"


template<typename T, typename I>
void submatrix_csx_csx<T,I>::set_matrix_src(MatrixCsxView_new * M_src_)
{
    M_src = M_src_;
}



template<typename T, typename I>
void submatrix_csx_csx<T,I>::set_bounds(size_t row_start_, size_t row_end_, size_t col_start_, size_t col_end_)
{
    if(M_src == nullptr) eslog::error("source matrix has not been set\n");

    row_start = row_start_;
    row_end = row_end_;
    col_start = col_start_;
    col_end = col_end_;
    num_rows = row_end - row_start;
    num_cols = col_end - col_start;

    if(row_start > row_end || row_end > M_src->nrows || col_start > col_end || col_end > M_src.ncols) eslog::error("wrong bounds\n");

    bounds_set = true;
}



template<typename T, typename I>
void submatrix_csx_csx<T,I>::setup()
{
    if(M_src == nullptr) eslog::error("source matrix has not been set\n");
    if(!bounds_set) eslog::error("bounds have not been set\n");

    size_t start_prim = 0;
    size_t end_prim = 0;
    size_t start_sec = 0;
    size_t end_sec = 0;
    if(M_src.order == 'R') {
        start_prim = start_row;
        end_prim = end_row;
        start_sec = start_col;
        end_sec = end_col;
    }
    if(M_src.order == 'C') {
        start_prim = start_col;
        end_prim = end_col;
        start_sec = start_row;
        end_sec = end_row;
    }

    nnz_output = 0;
    for(size_t ip = start_prim; ip < end_prim; ip++)
    {
        I start = M_src.ptrs[ip];
        I end = M_src.ptrs[ip+1];
        I i = start;
        while(i < end && M_src.idxs[i] < start_sec) {
            i++;
        }
        while(i < end && M_src.idxs[i] < end_sec) {
            nnz_output++;
            i++;
        }
    }

    setup_called = true;
}



template<typename T, typename I>
void submatrix_csx_csx<T,I>::get_output_matrix_nnz()
{
    if(!setup_called) eslog::error("setup has not been called\n");

    return nnz_output;
}



template<typename T, typename I>
void submatrix_csx_csx<T,I>::set_matrix_dst(MatrixCsxView_new * M_dst_)
{
    M_dst = M_dst_;
}



template<typename T, typename I>
void submatrix_csx_csx<T,I>::perform()
{
    if(M_src == nullptr) eslog::error("source matrix has not been set\n");
    if(M_dst == nullptr) eslog::error("destination matrix has not been set\n");
    if(!setup_called) eslog::error("setup has not been called\n");
    if(M_dst->nrows != num_rows || M_dst->ncols != num_cols) eslog::error("matrix sizes dont match\n");
    if(M_dst->order != M_src->order) eslog::error("matrix orders dont match\n");
    if(M_dst->nnz != nnz_output) eslog::error("wrong nnz in output matrix\n");

    size_t start_prim = 0;
    size_t end_prim = 0;
    size_t start_sec = 0;
    size_t end_sec = 0;
    if(M_src.order == 'R') {
        start_prim = start_row;
        end_prim = end_row;
        start_sec = start_col;
        end_sec = end_col;
    }
    if(M_src.order == 'C') {
        start_prim = start_col;
        end_prim = end_col;
        start_sec = start_row;
        end_sec = end_row;
    }

    size_t curr_nnz = 0;
    for(size_t ips = start_prim; ips < end_prim; ips++)
    {
        I ipd = ips - start_prim;
        M_dst.ptrs[ipd] = curr_nnz;
        I start = M_src.ptrs[ip];
        I end = M_src.ptrs[ip+1];
        I i = start;
        while(i < end) {
            I iss = M_src.idxs[i];
            if(iss >= start_sec) {
                break;
            }
            i++;
        }
        while(i < end) {
            I iss = M_src.idxs[i];
            if(iss >= end_sec) {
                break;
            }
            I isd = iss - start_sec;
            T val = M_src.vals[i];
            M_dst.idxs[curr_nnz] = isd;
            M_dst.vals[curr_nnz] = val;
            curr_nnz++;
            i++;
        }
    }
    M_dst.ptrs[curr_nnz] = nnz_output;
}



template<typename T, typename I>
void submatrix_csx_csx<T,I>::do_all(MatrixCsxView_new<T,I> * M_src, MatrixDenseView_new<T> * M_dst, size_t row_start, size_t row_end, size_t col_start, size_t col_end)
{
    submatrix_csx_csx<T,I> instance;
    instance.set_matrix_src(M_src);
    instance.set_matrix_dst(M_dst);
    instance.set_bounds(row_start, row_end, col_start, col_end);
    instance.setup();
    instance.perform();
}
