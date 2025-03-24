
#include "math/operations/submatrix_csx_csy_map.h"

#include "basis/utilities/stacktimer.h"



namespace espreso {
namespace math {
namespace operations {



template<typename T, typename I>
submatrix_csx_csy_map<T,I>::~submatrix_csx_csy_map()
{
    finalize();
}



template<typename T, typename I>
void submatrix_csx_csy_map<T,I>::set_matrix_src(MatrixCsxView_new<T,I> * M_src_)
{
    M_src = M_src_;
}



template<typename T, typename I>
void submatrix_csx_csy_map<T,I>::set_bounds(size_t row_start_, size_t row_end_, size_t col_start_, size_t col_end_)
{
    if(M_src == nullptr) eslog::error("source matrix has not been set\n");

    row_start = row_start_;
    row_end = row_end_;
    col_start = col_start_;
    col_end = col_end_;
    num_rows = row_end - row_start;
    num_cols = col_end - col_start;

    if(row_start > row_end || row_end > M_src->nrows || col_start > col_end || col_end > M_src->ncols) eslog::error("wrong bounds\n");

    set_bounds_called = true;
}



template<typename T, typename I>
void submatrix_csx_csy_map<T,I>::setup()
{
    stacktimer::push("submatrix_csx_csy_map::setup");

    if(M_src == nullptr) eslog::error("source matrix has not been set\n");
    if(!set_bounds_called) eslog::error("bounds have not been set\n");

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

    I * src_ptrs = M_src->ptrs;
    I * src_idxs = M_src->idxs;

    nnz_output = 0;
    for(I ip = start_prim; ip < end_prim; ip++)
    {
        I start = src_ptrs[ip];
        I end = src_ptrs[ip+1];
        I i = start;
        while(i < end && src_idxs[i] < start_sec) {
            i++;
        }
        while(i < end && src_idxs[i] < end_sec) {
            nnz_output++;
            i++;
        }
    }

    stacktimer::pop();

    setup_called = true;
}



template<typename T, typename I>
size_t submatrix_csx_csy_map<T,I>::get_output_matrix_nnz()
{
    if(!setup_called) eslog::error("setup has not been called\n");

    return nnz_output;
}



template<typename T, typename I>
void submatrix_csx_csy_map<T,I>::set_matrix_dst(MatrixCsxView_new<T,I> * M_dst_)
{
    M_dst = M_dst_;
}



template<typename T, typename I>
void submatrix_csx_csy_map<T,I>::perform_pattern()
{
    stacktimer::push("submatrix_csx_csy_map::perform_pattern");

    if(M_src == nullptr) eslog::error("source matrix has not been set\n");
    if(M_dst == nullptr) eslog::error("destination matrix has not been set\n");
    if(!setup_called) eslog::error("setup has not been called\n");
    if(M_dst->nrows != num_rows || M_dst->ncols != num_cols) eslog::error("matrix sizes dont match\n");
    if(M_dst->nnz != nnz_output) eslog::error("wrong nnz in output matrix\n");

    map.set(nnz_output, AllocatorCPU_new::get_singleton());
    map.alloc();

    if(M_src->order == M_dst->order) {
        perform_pattern_same_order();
    }
    else {
        perform_pattern_diff_order();
    }

    stacktimer::pop();

    perform_pattern_called = true;
}



template<typename T, typename I>
void submatrix_csx_csy_map<T,I>::perform_values()
{
    stacktimer::push("submatrix_csx_csy_map::perform_values");

    if(!perform_pattern_called) eslog::error("perform pattern has not been called\n");

    T * src_vals = M_src->vals;
    T * dst_vals = M_dst->vals;

    for(size_t i_dst = 0; i_dst < nnz_output; i_dst++) {
        I i_src = map.vals[i_dst];
        dst_vals[i_dst] = src_vals[i_src];
    }

    stacktimer::pop();
}



template<typename T, typename I>
void submatrix_csx_csy_map<T,I>::perform_all()
{
    perform_pattern();
    perform_values();
}



template<typename T, typename I>
void submatrix_csx_csy_map<T,I>::finalize()
{
    if(perform_pattern_called) {
        map.clear();
    }
    perform_pattern_called = false;
}



template<typename T, typename I>
void submatrix_csx_csy_map<T,I>::do_all(MatrixCsxView_new<T,I> * M_src, MatrixCsxView_new<T,I> * M_dst, size_t row_start, size_t row_end, size_t col_start, size_t col_end)
{
    submatrix_csx_csy_map<T,I> instance;
    instance.set_matrix_src(M_src);
    instance.set_matrix_dst(M_dst);
    instance.set_bounds(row_start, row_end, col_start, col_end);
    instance.setup();
    instance.perform_pattern();
    instance.perform_values();
    instance.finalize();
}



template<typename T, typename I>
void submatrix_csx_csy_map<T,I>::perform_pattern_same_order()
{
    if(M_src->order != M_dst->order) eslog::error("matrix orders dont match\n");
    
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

    I * src_ptrs = M_src->ptrs;
    I * src_idxs = M_src->idxs;
    I * dst_ptrs = M_dst->ptrs;
    I * dst_idxs = M_dst->idxs;

    size_t curr_nnz = 0;
    for(I ips = start_prim; ips < end_prim; ips++)
    {
        I ipd = ips - start_prim;
        dst_ptrs[ipd] = curr_nnz;
        I start = src_ptrs[ips];
        I end = src_ptrs[ips+1];
        I i = start;
        while(i < end) {
            I iss = src_idxs[i];
            if(iss >= start_sec) {
                break;
            }
            i++;
        }
        while(i < end) {
            I iss = src_idxs[i];
            if(iss >= end_sec) {
                break;
            }
            I isd = iss - start_sec;
            dst_idxs[curr_nnz] = isd;
            map.vals[curr_nnz] = i;
            curr_nnz++;
            i++;
        }
    }
    dst_ptrs[M_dst->get_size_primary()] = nnz_output;
}



template<typename T, typename I>
void submatrix_csx_csy_map<T,I>::perform_pattern_diff_order()
{
    if(M_src->order == M_dst->order) eslog::error("matrix orders must not match\n");

    // actually more similar to convert_csx_csy than submatrix
    // ipd = index primary destination, ...
    // I src_size_primary = M_src->get_size_primary();
    // size_t src_size_secdary = M_src->get_size_secdary();
    I dst_size_primary = M_dst->get_size_primary();
    // size_t dst_size_secdary = M_dst->get_size_secdary();
    // size_t nnz = M_src->nnz;
    I * src_ptrs = M_src->ptrs;
    I * src_idxs = M_src->idxs;
    I * dst_ptrs = M_dst->ptrs;
    I * dst_idxs = M_dst->idxs;

    I src_primary_start = ((M_src->order == 'R') ? row_start : col_start);
    I src_primary_end = ((M_src->order == 'R') ? row_end : col_end);
    I src_secdary_start = ((M_src->order == 'R') ? col_start : row_start);
    I src_secdary_end = ((M_src->order == 'R') ? col_end : row_end);

    // initialize nnz per dst primary to 0
    for(I ipd = 0; ipd <= dst_size_primary; ipd++) {
        dst_ptrs[ipd] = 0;
    }

    // calculate dst nnz per primary
    for(I ips = src_primary_start; ips < src_primary_end; ips++) {
        I start = src_ptrs[ips];
        I end = src_ptrs[ips+1];
        I i = start;
        I iss;
        while(i < end && (iss = src_idxs[i]) < src_secdary_start) {
            i++;
        }
        while(i < end && (iss = src_idxs[i]) < src_secdary_end) {
            I ipd = iss - src_secdary_start;
            dst_ptrs[ipd]++;
            i++;
        }
    }

    // exclusive cumulative sum
    I curr = 0;
    for(I ipd = 0; ipd <= dst_size_primary; ipd++)
    {
        I tmp = dst_ptrs[ipd];
        dst_ptrs[ipd] = curr;
        curr += tmp;
    }

    // fill dst idxs and map
    for(I ips = src_primary_start; ips < src_primary_end; ips++)
    {
        I start = src_ptrs[ips];
        I end = src_ptrs[ips+1];
        I i_src = start;
        I iss;
        while(i_src < end && (iss = src_idxs[i_src]) < src_secdary_start) {
            i_src++;
        }
        while(i_src < end && (iss = src_idxs[i_src]) < src_secdary_end) {
            I ipd = iss - src_secdary_start;
            I i_dst = dst_ptrs[ipd];
            dst_ptrs[ipd]++;
            I isd = ips - src_primary_start;
            dst_idxs[i_dst] = isd;
            map.vals[i_dst] = i_src;
            i_src++;
        }
    }

    // fix (shift) dst ptrs
    curr = 0;
    for(I ipd = 0; ipd <= dst_size_primary; ipd++)
    {
        I tmp = dst_ptrs[ipd];
        dst_ptrs[ipd] = curr;
        curr = tmp;
    }
}



#define INSTANTIATE_T_I(T,I) \
template class submatrix_csx_csy_map<T,I>;

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
