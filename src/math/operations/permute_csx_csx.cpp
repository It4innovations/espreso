
#include "math/operations/permute_csx_csx.h"

#include "math/primitives_new/allocator_new.h"
#include "math/operations/copy_csx.h"
#include "math/primitives_new/vector_dense_data_new.h"
#include "basis/utilities/stacktimer.h"



namespace espreso {
namespace math {
namespace operations {



template<typename T, typename I>
void permute_csx_csx<T,I>::set_matrix_src(MatrixCsxView_new<T,I> * M_src_)
{
    if(M_src != nullptr) eslog::error("matrix M_src is already set\n");

    M_src = M_src_;
}



template<typename T, typename I>
void permute_csx_csx<T,I>::set_matrix_dst(MatrixCsxView_new<T,I> * M_dst_)
{
    if(M_dst != nullptr) eslog::error("matrix M_dst is already set\n");

    M_dst = M_dst_;
}



template<typename T, typename I>
void permute_csx_csx<T,I>::set_perm_rows(PermutationView_new<I> * perm_rows_)
{
    if(perm_rows != nullptr) eslog::error("perm_rows is already set\n");

    perm_rows = perm_rows_;
}



template<typename T, typename I>
void permute_csx_csx<T,I>::set_perm_cols(PermutationView_new<I> * perm_cols_)
{
    if(perm_cols != nullptr) eslog::error("perm_cols is already set\n");

    perm_cols = perm_cols_;
}



template<typename T, typename I>
void permute_csx_csx<T,I>::perform()
{
    stacktimer::push("permute_csx_csx::perform");

    if(M_src == nullptr) eslog::error("source matrix is not set\n");
    if(M_dst == nullptr) eslog::error("destination matrix is not set\n");
    if(!M_src->ator->is_data_accessible_cpu()) eslog::error("source matrix must be cpu-accessible\n");
    if(!M_dst->ator->is_data_accessible_cpu()) eslog::error("destination matrix must be cpu-accessible\n");
    if(M_src->order != M_dst->order) eslog::error("matrix orders dont match\n");
    if(M_src->nrows != M_dst->nrows || M_src->ncols != M_dst->ncols || M_src->nnz != M_dst->nnz) eslog::error("matrix sizes dont match\n");
    if(perm_rows != nullptr && perm_rows->size != M_src->nrows) eslog::error("wrong row perm size\n");
    if(perm_cols != nullptr && perm_cols->size != M_src->ncols) eslog::error("wrong col perm size\n");

    if(perm_rows == nullptr && perm_cols == nullptr) {
        copy_csx<T,I>::do_all(M_src, M_dst);
    }
    if(perm_rows != nullptr && perm_cols == nullptr) {
        if(M_src->order == 'R') {
            perform_primary(*perm_rows);
        }
        if(M_src->order == 'C') {
            perform_secdary(*perm_rows);
        }
    }
    if(perm_rows == nullptr && perm_cols != nullptr) {
        if(M_src->order == 'R') {
            perform_secdary(*perm_cols);
        }
        if(M_src->order == 'C') {
            perform_primary(*perm_cols);
        }
    }
    if(perm_rows != nullptr && perm_cols != nullptr) {
        if(M_src->order == 'R') {
            perform_both(*perm_rows, *perm_cols);
        }
        if(M_src->order == 'C') {
            perform_both(*perm_cols, *perm_rows);
        }
    }

    stacktimer::pop();
}



template<typename T, typename I>
void permute_csx_csx<T,I>::do_all(MatrixCsxView_new<T,I> * M_src, MatrixCsxView_new<T,I> * M_dst, PermutationView_new<I> * perm_rows, PermutationView_new<I> * perm_cols)
{
    permute_csx_csx<T,I> instance;
    instance.set_matrix_src(M_src);
    instance.set_matrix_dst(M_dst);
    instance.set_perm_rows(perm_rows);
    instance.set_perm_cols(perm_cols);
    instance.perform();
}



template<typename T, typename I>
void permute_csx_csx<T,I>::perform_primary(PermutationView_new<I> & perm)
{
    I * src_ptrs = M_src->ptrs;
    I * src_idxs = M_src->idxs;
    T * src_vals = M_src->vals;
    I * dst_ptrs = M_dst->ptrs;
    I * dst_idxs = M_dst->idxs;
    T * dst_vals = M_dst->vals;
    size_t size_primary = M_src->get_size_primary();
    // size_t size_secdary = M_src->get_size_secdary();

    I i_dst = 0;
    for(size_t ipd = 0; ipd < size_primary; ipd++) {
        dst_ptrs[ipd] = i_dst;
        I ips = perm.dst_to_src[ipd];
        I start = src_ptrs[ips];
        I end = src_ptrs[ips+1];
        for(I i_src = start; i_src < end; i_src++) {
            I is = src_idxs[i_src];
            T val = src_vals[i_src];
            dst_idxs[i_dst] = is;
            dst_vals[i_dst] = val;
            i_dst++;
        }
    }
    dst_ptrs[size_primary] = i_dst;
}



template<typename T, typename I>
void permute_csx_csx<T,I>::perform_secdary(PermutationView_new<I> & perm)
{
    I * src_ptrs = M_src->ptrs;
    I * src_idxs = M_src->idxs;
    T * src_vals = M_src->vals;
    I * dst_ptrs = M_dst->ptrs;
    I * dst_idxs = M_dst->idxs;
    T * dst_vals = M_dst->vals;
    size_t size_primary = M_src->get_size_primary();
    // size_t size_secdary = M_src->get_size_secdary();

    std::copy_n(src_ptrs, size_primary + 1, dst_ptrs);

    I longest_secdary = 0;
    for(size_t ip = 0; ip < size_primary; ip++) {
        I start = src_ptrs[ip];
        I end = src_ptrs[ip+1];
        I length = end - start;
        longest_secdary = std::max(longest_secdary, length);
    }

    struct isd_val { I isd; T val; };

    VectorDenseData_new<isd_val> ivs;
    ivs.set(longest_secdary, AllocatorCPU_new::get_singleton());
    ivs.alloc();

    for(size_t ip = 0; ip < size_primary; ip++) {
        I start = src_ptrs[ip];
        I end = src_ptrs[ip+1];
        for(I i = start; i < end; i++) {
            I iss = src_idxs[i];
            I isd = perm.src_to_dst[iss];
            T val = src_vals[i];
            ivs.vals[i - start] = isd_val{isd, val};
        }
        std::sort(ivs.vals, ivs.vals + end - start, [](const isd_val & l, const isd_val & r){return l.isd < r.isd;});
        for(I i = start; i < end; i++) {
            I isd = ivs.vals[i - start].isd;
            T val = ivs.vals[i - start].val;
            dst_idxs[i] = isd;
            dst_vals[i] = val;
        }
    }

    ivs.clear();
}



template<typename T, typename I>
void permute_csx_csx<T,I>::perform_both(PermutationView_new<I> & perm_primary, PermutationView_new<I> & perm_secdary)
{
    I * src_ptrs = M_src->ptrs;
    I * src_idxs = M_src->idxs;
    T * src_vals = M_src->vals;
    I * dst_ptrs = M_dst->ptrs;
    I * dst_idxs = M_dst->idxs;
    T * dst_vals = M_dst->vals;
    size_t size_primary = M_src->get_size_primary();
    // size_t size_secdary = M_src->get_size_secdary();

    I longest_secdary = 0;
    for(size_t ip = 0; ip < size_primary; ip++) {
        I start = src_ptrs[ip];
        I end = src_ptrs[ip+1];
        I length = end - start;
        longest_secdary = std::max(longest_secdary, length);
    }

    struct isd_val { I isd; T val; };

    VectorDenseData_new<isd_val> ivs;
    ivs.set(longest_secdary, AllocatorCPU_new::get_singleton());
    ivs.alloc();

    I i_dst = 0;
    for(size_t ipd = 0; ipd < size_primary; ipd++) {
        dst_ptrs[ipd] = i_dst;
        I ips = perm_primary.dst_to_src[ipd];
        I start = src_ptrs[ips];
        I end = src_ptrs[ips+1];
        for(I i_src = start; i_src < end; i_src++) {
            I iss = src_idxs[i_src];
            I isd = perm_secdary.src_to_dst[iss];
            T val = src_vals[i_src];
            ivs.vals[i_src - start] = isd_val{isd, val};
        }
        std::sort(ivs.vals, ivs.vals + end - start, [](const isd_val & l, const isd_val & r){return l.isd < r.isd;});
        for(I i_src = start; i_src < end; i_src++) {
            I isd = ivs.vals[i_src - start].isd;
            T val = ivs.vals[i_src - start].val;
            dst_idxs[i_dst] = isd;
            dst_vals[i_dst] = val;
            i_dst++;
        }
    }
    dst_ptrs[size_primary] = i_dst;

    ivs.clear();
}



#define INSTANTIATE_T_I(T,I) \
template class permute_csx_csx<T,I>;

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
