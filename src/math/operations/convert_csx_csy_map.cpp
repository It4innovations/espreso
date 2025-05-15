
#include "math/operations/convert_csx_csy_map.h"

#include "math/operations/copy_csx.h"
#include "math/primitives_new/allocator_new.h"
#include "basis/utilities/stacktimer.h"



namespace espreso {
namespace math {
namespace operations {



template<typename T, typename I>
void convert_csx_csy_map<T,I>::set_matrix_src(MatrixCsxView_new<T,I> * M_src_)
{
    if(M_src != nullptr) eslog::error("matrix M_src is already set\n");

    M_src = M_src_;
}



template<typename T, typename I>
void convert_csx_csy_map<T,I>::set_matrix_dst(MatrixCsxView_new<T,I> * M_dst_)
{
    if(M_dst != nullptr) eslog::error("matrix M_dst is already set\n");

    M_dst = M_dst_;
}



template<typename T, typename I>
void convert_csx_csy_map<T,I>::perform_pattern()
{
    stacktimer::push("convert_csx_csy_map::perform_pattern");

    if(perform_pattern_called) eslog::error("pattern computation was already performed\n");
    if(M_src == nullptr) eslog::error("source matrix is not set\n");
    if(M_dst == nullptr) eslog::error("destination matrix is not set\n");
    if(!M_src->ator->is_data_accessible_cpu()) eslog::error("source matrix must be cpu-accessible\n");
    if(!M_dst->ator->is_data_accessible_cpu()) eslog::error("destination matrix must be cpu-accessible\n");
    if(M_dst->nrows != M_src->nrows || M_dst->ncols != M_src->ncols || M_dst->nnz != M_src->nnz) eslog::error("matrix sizes dont match\n");

    if(M_src->order == M_dst->order) {
        std::copy_n(M_src->ptrs, M_src->get_size_primary() + 1, M_dst->ptrs);
        std::copy_n(M_src->idxs, M_src->nnz, M_dst->idxs);
        stacktimer::pop();
        perform_pattern_called = true;
        return;
    }

    // use terminology for CSR->CSC, the other way it works equally

    size_t nrows = M_src->get_size_primary();
    size_t ncols = M_src->get_size_secdary();
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
    for(size_t c = 0; c <= ncols; c++) {
        dst_ptrs[c] = 0;
    }

    // calculate dst nnz per col
    for(size_t i = 0; i < nnz; i++) {
        I c = src_idxs[i];
        dst_ptrs[c]++;
    }

    // exclusive cumulative sum
    I curr = 0;
    for(size_t c = 0; c <= ncols; c++)
    {
        I tmp = dst_ptrs[c];
        dst_ptrs[c] = curr;
        curr += tmp;
    }

    // fill dst idxs and map
    for(size_t r = 0; r < nrows; r++)
    {
        I start = src_ptrs[r];
        I end = src_ptrs[r+1];
        for(I i_src = start; i_src < end; i_src++)
        {
            I c = src_idxs[i_src];
            I i_dst = dst_ptrs[c];
            dst_ptrs[c]++;
            dst_idxs[i_dst] = r;
            map_vals[i_dst] = i_src;
        }
    }

    // fix (shift) dst ptrs
    curr = 0;
    for(size_t c = 0; c <= ncols; c++)
    {
        I tmp = dst_ptrs[c];
        dst_ptrs[c] = curr;
        curr = tmp;
    }

    stacktimer::pop();

    perform_pattern_called = true;
}



template<typename T, typename I>
void convert_csx_csy_map<T,I>::perform_values()
{
    stacktimer::push("convert_csx_csy_map::perform_values");

    if(!perform_pattern_called) eslog::error("pattern computation has not been performed\n");
    if(M_src == nullptr) eslog::error("source matrix is not set\n");
    if(M_dst == nullptr) eslog::error("destination matrix is not set\n");
    if(M_dst->nrows != M_src->nrows || M_dst->ncols != M_src->ncols || M_dst->nnz != M_src->nnz) eslog::error("matrix sizes dont match\n");

    if(M_src->order == M_dst->order) {
        std::copy_n(M_src->vals, M_src->nnz, M_dst->vals);
        stacktimer::pop();
        return;
    }

    T * src_vals = M_src->vals;
    T * dst_vals = M_dst->vals;
    I * map_vals = map.vals;
    I nnz = M_src->nnz;

    for(I i_dst = 0; i_dst < nnz; i_dst++) {
        I i_src = map_vals[i_dst];
        dst_vals[i_dst] = src_vals[i_src];
    }

    stacktimer::pop();
}



template<typename T, typename I>
void convert_csx_csy_map<T,I>::perform_all()
{
    perform_pattern();
    perform_values();
}



template<typename T, typename I>
void convert_csx_csy_map<T,I>::do_all(MatrixCsxView_new<T,I> * M_src, MatrixCsxView_new<T,I> * M_dst)
{
    convert_csx_csy_map<T,I> instance;
    instance.set_matrix_src(M_src);
    instance.set_matrix_dst(M_dst);
    instance.perform_pattern();
    instance.perform_values();
}



#define INSTANTIATE_T_I(T,I) \
template class convert_csx_csy_map<T,I>;

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
