
#include "math/operations/submatrix_csx_dny.h"

#include "basis/utilities/stacktimer.h"



namespace espreso {
namespace math {
namespace operations {



template<typename T, typename I>
void submatrix_csx_dny<T,I>::set_matrix_src(MatrixCsxView_new<T,I> * M_src_)
{
    M_src = M_src_;
}



template<typename T, typename I>
void submatrix_csx_dny<T,I>::set_matrix_dst(MatrixDenseView_new<T> * M_dst_)
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
    num_rows = row_end - row_start;
    num_cols = col_end - col_start;

    if(row_start > row_end || row_end > M_src->nrows || col_start > col_end || col_end > M_src->ncols) eslog::error("wrong bounds\n");

    bound_set = true;
}



template<typename T, typename I>
void submatrix_csx_dny<T,I>::perform_zerofill()
{
    if(M_src == nullptr) eslog::error("source matrix is not set\n");
    if(M_dst == nullptr) eslog::error("destination matrix is not set\n");
    if(!bound_set) eslog::error("bounds are not set\n");
    if(M_dst->nrows != num_rows || M_dst->ncols != num_cols) eslog::error("wrong output matrix size\n");

    stacktimer::push("submatrix_csx_dny::perform_zerofill");
    
    size_t size_primary = M_dst->get_size_primary();
    size_t size_secdary = M_dst->get_size_secdary();
    for(size_t i = 0; i < size_primary; i++) {
        std::fill_n(M_dst->vals + i * M_dst->ld, size_secdary, T{0});
    }

    stacktimer::pop();

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

    stacktimer::push("submatrix_csx_dny::perform_copyvals");

    I start_prim = 0;
    I end_prim = 0;
    I start_sec = 0;
    I end_sec = 0;
    if(M_src->order == 'R') {
        start_prim = row_start;
        end_prim = row_end;
        start_sec = col_start;
        end_sec = col_end;
    }
    if(M_src->order == 'C') {
        start_prim = col_start;
        end_prim = col_end;
        start_sec = row_start;
        end_sec = row_end;
    }

    size_t stride_prim = ((M_src->order == M_dst->order) ? M_dst->ld : 1);
    size_t stride_sec = ((M_src->order == M_dst->order) ? 1 : M_dst->ld);

    I * srcptrs = M_src->ptrs;
    I * srcidxs = M_src->idxs;
    T * srcvals = M_src->vals;
    T * dstvals = M_dst->vals;

    for(I ips = start_prim; ips < end_prim; ips++)
    {
        size_t ipd = ips - start_prim;
        I start = srcptrs[ips];
        I end = srcptrs[ips+1];
        I i = start;
        while(i < end && srcidxs[i] < start_sec) {
            i++;
        }
        while(i < end && srcidxs[i] < end_sec) {
            I iss = srcidxs[i];
            I isd = iss - start_sec;
            T v = srcvals[i];
            dstvals[ipd * stride_prim + isd * stride_sec] = v;
            i++;
        }
    }

    stacktimer::pop();
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
        INSTANTIATE_T(std::complex<double>)

            INSTANTIATE

        #undef INSTANTIATE
    #undef INSTANTIATE_T
#undef INSTANTIATE_T_I



}
}
}
