
#include "math/operations/permute_dnx_dnx.h"

#include "math/primitives_new/allocator_new.h"
#include "math/primitives_new/matrix_dense_data_new.h"
#include "math/operations/complete_dnx_dnx.h"
#include "math/operations/copy_dnx.h"
#include "basis/utilities/stacktimer.h"



namespace espreso {
namespace math {
namespace operations {



template<typename T, typename I>
void permute_dnx_dnx<T,I>::set_matrix_src(MatrixDenseView_new<T> * M_src_)
{
    if(M_src != nullptr) eslog::error("matrix M_src is already set\n");

    M_src = M_src_;
}



template<typename T, typename I>
void permute_dnx_dnx<T,I>::set_matrix_dst(MatrixDenseView_new<T> * M_dst_)
{
    if(M_dst != nullptr) eslog::error("matrix M_dst is already set\n");

    M_dst = M_dst_;
}



template<typename T, typename I>
void permute_dnx_dnx<T,I>::set_perm_vector_rows(PermutationView_new<I> * perm_rows_)
{
    if(perm_rows != nullptr) eslog::error("perm_rows is already set\n");

    perm_rows = perm_rows_;
}



template<typename T, typename I>
void permute_dnx_dnx<T,I>::set_perm_vector_cols(PermutationView_new<I> * perm_cols_)
{
    if(perm_cols != nullptr) eslog::error("perm_cols is already set\n");

    perm_cols = perm_cols_;
}



template<typename T, typename I>
void permute_dnx_dnx<T,I>::perform()
{
    stacktimer::push("permute_dnx_dnx::perform");

    if(M_src == nullptr) eslog::error("source matrix is not set\n");
    if(M_dst == nullptr) eslog::error("destination matrix is not set\n");
    if(!M_src->ator->is_data_accessible_cpu()) eslog::error("source matrix must be cpu-accessible\n");
    if(!M_dst->ator->is_data_accessible_cpu()) eslog::error("destination matrix must be cpu-accessible\n");
    if(M_src->order != M_dst->order) eslog::error("matrix orders dont match\n");
    if(!is_uplo_equal(M_src->prop.uplo, M_dst->prop.uplo)) eslog::error("matrix uplo does not match\n");
    if(M_src->nrows != M_dst->nrows || M_src->ncols != M_dst->ncols) eslog::error("matrix sizes dont match\n");
    if(perm_rows != nullptr && perm_rows->size != M_src->nrows) eslog::error("wrong row perm size\n");
    if(perm_cols != nullptr && perm_cols->size != M_src->ncols) eslog::error("wrong col perm size\n");

    if(M_dst->prop.uplo == 'U' || M_dst->prop.uplo == 'L') {
        if(M_src->nrows != M_src->ncols) eslog::error("uplo can be set only for square matrices\n");

        MatrixDenseData_new<T> tmp1;
        tmp1.set(M_src->nrows, M_src->ncols, M_src->order, AllocatorCPU_new::get_singleton());
        tmp1.alloc();

        MatrixDenseData_new<T> tmp2;
        tmp2.set(M_src->nrows, M_src->ncols, M_src->order, AllocatorCPU_new::get_singleton());
        tmp2.alloc();

        complete_dnx_dnx<T>::do_all(M_src, &tmp1, is_hermitian<T>(M_src->prop.symm));
        permute_dnx_dnx<T,I>::do_all(&tmp1, &tmp2, perm_rows, perm_cols);
        tmp2.prop.uplo = M_dst->prop.uplo;
        copy_dnx<T>::do_all(&tmp2, M_dst, false);
    }
    else {
        if(M_src == M_dst) {
            MatrixDenseData_new<T> tmp;
            tmp.set(M_src->nrows, M_src->ncols, M_src->order, AllocatorCPU_new::get_singleton());
            tmp.alloc();

            permute_dnx_dnx<T,I>::do_all(M_src, &tmp, perm_rows, perm_cols);
            copy_dnx<T>::do_all(&tmp, M_dst, false);
        }
        else {
            if(perm_rows == nullptr && perm_cols == nullptr) {
                copy_dnx<T>::do_all(M_src, M_dst);
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
    }

    stacktimer::pop();
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
        I ips = perm_primary.dst_to_src[ipd];
        T * src_sub = src_vals + ips * ld_src;
        T * dst_sub = dst_vals + ipd * ld_dst;
        for(size_t isd = 0; isd < size_secdary; isd++) {
            I iss = perm_secdary.dst_to_src[isd];
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
        INSTANTIATE_T(std::complex<double>)

            INSTANTIATE

        #undef INSTANTIATE
    #undef INSTANTIATE_T
#undef INSTANTIATE_T_I



}
}
}
