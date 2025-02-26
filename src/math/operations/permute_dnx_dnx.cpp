
#include "math/operations/permute_dnx_dnx.h"



namespace espreso {
namespace math {
namespace operations {



template<typename T, typename I>
void permute_dnx_dnx<T,I>::set_matrix_src(MatrixDenseView_new<T> * M_src_)
{
    M_src = M_src_;
}



template<typename T, typename I>
void permute_dnx_dnx<T,I>::set_matrix_dst(MatrixDenseView_new<T> * M_dst_)
{
    M_dst = M_dst_;
}



template<typename T, typename I>
void permute_dnx_dnx<T,I>::set_perm_vector_rows(PermutationView_new<I> * perm_rows_)
{
    perm_rows = perm_rows_;
}



template<typename T, typename I>
void permute_dnx_dnx<T,I>::set_perm_vector_cols(PermutationView_new<I> * perm_cols_)
{
    perm_cols = perm_cols_;
}



template<typename T, typename I>
void permute_dnx_dnx<T,I>::perform()
{
    if(M_src == nullptr) eslog::error("source matrix is not set\n");
    if(M_dst == nullptr) eslog::error("destination matrix is not set\n");
    if(M_src->order != M_dst->order) eslog::error("matrix orders dont match\n");
    if(M_src->nrows != M_dst->nrows || M_src->ncols != M_dst->ncols) eslog::error("matrix sizes dont match\n");
    if(perm_rows != nullptr && perm_rows->size != M_src->nrows) eslog::error("wrong row perm size\n");
    if(perm_cols != nullptr && perm_cols->size != M_src->ncols) eslog::error("wrong col perm size\n");

    if(perm_rows == nullptr && perm_cols == nullptr) {
        copy_dense<T>::do_all(M_src, M_dst);
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
}



template<typename T, typename I>
void permute_dnx_dnx<T,I>::do_all(MatrixDenseView_new<T> * M_src, MatrixDenseView_new<T> * M_dst, PermutationView_new<I> * perm_rows, PermutationView_new<I> * perm_cols)
{
    permute_dnx_dnx<T,I> instance;
    instance.set_matrix_src(M_src);
    instance.set_matrix_dst(M_dst);
    instance.set_perm_vector_rows(perm_rows);
    instance.set_perm_vector_cols(perm_cols);
    instance.perform();
}



template<typename T, typename I>
void permute_dnx_dnx<T,I>::perform_primary(PermutationView_new<I> & perm)
{
    T * src_vals = M_src->vals;
    T * dst_vals = M_dst->vals;
    size_t size_primary = M_dst->get_size_primary();
    size_t size_secdary = M_dst->get_size_secdary();
    size_t ld_src = M_src->ld;
    size_t ld_dst = M_dst->ld;

    for(size_t ipd = 0; ipd < size_primary; ipd++) {
        I ips = perm.dst_to_src[ipd];
        std::copy_n(src_vals + ips * ld_src, size_secdary, dst_vals + ipd * ld_dst);
    }
}



template<typename T, typename I>
void permute_dnx_dnx<T,I>::perform_secdary(PermutationView_new<I> & perm)
{
    T * src_vals = M_src->vals;
    T * dst_vals = M_dst->vals;
    size_t size_primary = M_dst->get_size_primary();
    size_t size_secdary = M_dst->get_size_secdary();
    size_t ld_src = M_src->ld;
    size_t ld_dst = M_dst->ld;

    for(size_t ip = 0; ip < size_primary; ip++) {
        T * src_sub = src_vals + ip * ld_src;
        T * dst_sub = dst_vals + ip * ld_dst;
        for(size_t isd = 0; isd < size_secdary; isd++) {
            I iss = perm.dst_to_src[isd];
            dst_sub[isd] = src_sub[iss];
        }
    }
}



template<typename T, typename I>
void permute_dnx_dnx<T,I>::perform_both(PermutationView_new<I> & perm_primary, PermutationView_new<I> & perm_secdary)
{
    T * src_vals = M_src->vals;
    T * dst_vals = M_dst->vals;
    size_t size_primary = M_dst->get_size_primary();
    size_t size_secdary = M_dst->get_size_secdary();
    size_t ld_src = M_src->ld;
    size_t ld_dst = M_dst->ld;

    for(size_t ipd = 0; ipd < size_primary; ipd++) {
        I ips = perm.dst_to_src[ipd];
        T * src_sub = src_vals + ips * ld_src;
        T * dst_sub = dst_vals + ipd * ld_dst;
        for(size_t isd = 0; isd < size_secdary; isd++) {
            I iss = perm.dst_to_src[isd];
            dst_sub[isd] = src_sub[iss];
        }        
    }
}



#define INSTANTIATE_T_I(T,I) \
template class permute_dnx_dnx<T,I>;

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
