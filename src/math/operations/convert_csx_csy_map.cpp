
#include "math/operations/convert_csx_csy_map.h"

#include "math/operations/copy_csx.h"



template<typename T, typename I>
convert_csx_csy_map<T,I>::~convert_csx_csy_map()
{
    finalize();
}



template<typename T, typename I>
void convert_csx_csy_map<T,I>::set_matrix_src(MatrixCsxView_new<T,I> * M_src_)
{
    M_src = M_src_;
}



template<typename T, typename I>
void convert_csx_csy_map<T,I>::set_matrix_dst(MatrixCsxView_new<T,I> * M_dst_)
{
    M_dst = M_dst_;
}



template<typename T, typename I>
void convert_csx_csy_map<T,I>::perform_pattern()
{
    if(perform_pattern_called) eslog::error("pattern computation was already performed\n");
    if(M_src == nullptr) eslog::error("source matrix is not set\n");
    if(M_dst == nullptr) eslog::error("destination matrix is not set\n");
    if(M_dst->nrows != M_src->nrows || M_dst->ncols != M_src->ncols || M_dst->nnz != M_src->nnz) eslog::error("matrix sizes dont match\n");

    if(M_src->order == M_dst->order) {
        std::copy_n(M_src->ptrs, M_src->get_primary_size() + 1, M_dst->ptrs);
        std::copy_n(M_src->idxs, M_src->nnz, M_dst->idxs);
        return;
    }

    // use terminology for CSR->CSC, the other way it works equally

    size_t nrows_src = M_src->get_primary_size();
    size_t ncols_src = M_src->get_secondary_size();
    size_t nrows_dst = M_dst->get_primary_size();
    size_t ncols_dst = M_dst->get_secondary_size();
    size_t nnz = M_src->nnz;
    I * src_ptrs = M_src->ptrs;
    I * src_idxs = M_src->idxs;
    I * dst_ptrs = M_dst->ptrs;
    I * dst_idxs = M_dst->idxs;

    // initialize map
    map.set(nnz, AllocatorCPU_new::get_singleton());
    map.alloc();
    I * map_vals = map.vals;
    
    // initialize nnz per dst col to 0
    for(size_t c = 0; c <= ncols_dst; c++) {
        dst_ptrs[c] = 0;
    }

    // calculate dst nnz per col
    for(size_t i = 0; i < nnz; i++) {
        I c = src_idxs[i];
        dst_ptrs[c]++;
    }

    // exclusive cumulative sum
    I curr = 0;
    for(size_t c = 0; c <= ncols_dst; c++)
    {
        I tmp = dst_ptrs[c];
        dst_ptrs[c] = curr;
        curr += tmp;
    }

    // fill dst idxs and map
    for(size_t r = 0; r < nrows_src; r++)
    {
        I start = src_ptrs[r_in];
        I end = src_ptrs[r_in+1];
        for(I i_src = start; i_src < end; i_src++)
        {
            I c = src_idxs[i];
            I i_dst = M_dst.ptrs[c];
            dst_ptrs[c]++;
            dst_idxs[i_dst] = r;
            map_vals[i_dst] = i_src;
        }
    }

    // fix (shift) dst ptrs
    curr = 0;
    for(size_t c = 0; c <= ncols_dst; c++)
    {
        I tmp = dst_ptrs[c];
        dst_ptrs[c] = curr;
        curr += tmp;
    }

    perform_pattern_called = true;
}



template<typename T, typename I>
void convert_csx_csy_map<T,I>::perform_values()
{
    if(perform_pattern_called) eslog::error("pattern computation has not been performed\n");
    if(M_src == nullptr) eslog::error("source matrix is not set\n");
    if(M_dst == nullptr) eslog::error("destination matrix is not set\n");
    if(M_dst->nrows != M_src->nrows || M_dst->ncols != M_src->ncols || M_dst->nnz != M_src->nnz) eslog::error("matrix sizes dont match\n");

    if(M_src->order == M_dst->order) {
        std::copy_n(M_src->vals, M_src->nnz, M_dst->vals);
        return;
    }

    I * src_vals = M_src->vals;
    I * dst_vals = M_dst->vals;
    I * map_vals = map.vals;
    I nnz = M_src->nnz;

    for(I i_dst = 0; i_dst < nnz; i_dst++) {
        I i_src = map_vals[i_dst];
        dst_vals[i_dst] = src_vals[i_src];
    }
}



template<typename T, typename I>
void convert_csx_csy_map<T,I>::perform_all()
{
    perform_pattern();
    perform_values();
}



template<typename T, typename I>
void convert_csx_csy_map<T,I>::finalize()
{
    map.clear();
}



template<typename T, typename I>
void convert_csx_csy_map<T,I>::do_all(MatrixCsxView_new<T,I> * M_src, MatrixCsxView_new<T,I> * M_dst)
{
    convert_csx_csy_map<T,I> instance;
    instance.set_matrix_src(M_src);
    instance.set_matrix_dst(M_dst);
    instance.perform_pattern();
    instance.perform_values();
    instance.finalize();
}
