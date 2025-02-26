
#include "math/operations/convert_csx_csy.h"

#include "math/operations/copy_csx.h"



namespace espreso {
namespace math {
namespace operations {



template<typename T, typename I>
void convert_csx_csy<T,I>::set_matrix_src(MatrixCsxView_new<T,I> * M_src_)
{
    M_src = M_src_;
}



template<typename T, typename I>
void convert_csx_csy<T,I>::set_matrix_dst(MatrixCsxView_new<T,I> * M_dst_)
{
    M_dst = M_dst_;
}



template<typename T, typename I>
void convert_csx_csy<T,I>::perform()
{
    if(M_src == nullptr) eslog::error("source matrix is not set\n");
    if(M_dst == nullptr) eslog::error("destination matrix is not set\n");
    if(M_dst->nrows != M_src->nrows || M_dst->ncols != M_src->ncols || M_dst->nnz != M_src->nnz) eslog::error("matrix sizes dont match\n");

    if(M_src->order == M_dst->order) {
        copy_csx<T,I>::do_all(M_src, M_dst);
        return;
    }

    // use terminology for CSR->CSC, the other way it works equally

    size_t nrows_src = M_src->get_size_primary();
    // size_t ncols_src = M_src->get_size_secdary();
    // size_t nrows_dst = M_dst->get_size_secdary();
    size_t ncols_dst = M_dst->get_size_primary();
    size_t nnz = M_src->nnz;
    I * src_ptrs = M_src->ptrs;
    I * src_idxs = M_src->idxs;
    T * src_vals = M_src->vals;
    I * dst_ptrs = M_dst->ptrs;
    I * dst_idxs = M_dst->idxs;
    T * dst_vals = M_dst->vals;
    
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

    // fill dst idxs and vals
    for(size_t r = 0; r < nrows_src; r++)
    {
        I start = src_ptrs[r];
        I end = src_ptrs[r+1];
        for(I i_src = start; i_src < end; i_src++)
        {
            I c = src_idxs[i_src];
            T v = src_vals[i_src];
            I i_dst = dst_ptrs[c];
            dst_ptrs[c]++;
            dst_idxs[i_dst] = r;
            dst_vals[i_dst] = v;
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
}



template<typename T, typename I>
void convert_csx_csy<T,I>::do_all(MatrixCsxView_new<T,I> * M_src, MatrixCsxView_new<T,I> * M_dst)
{
    convert_csx_csy<T,I> instance;
    instance.set_matrix_src(M_src);
    instance.set_matrix_dst(M_dst);
    instance.perform();
}



#define INSTANTIATE_T_I(T,I) \
template class convert_csx_csy<T,I>;

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
