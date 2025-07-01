
#include "math/operations/permute_csx_dny.h"

#include "math/operations/fill_dnx.h"
#include "basis/utilities/stacktimer.h"



namespace espreso {
namespace math {
namespace operations {



template<typename T, typename I>
void permute_csx_dny<T,I>::set_matrix_src(MatrixCsxView_new<T,I> * M_src_)
{
    if(M_src != nullptr) eslog::error("matrix M_src is already set\n");

    M_src = M_src_;
}



template<typename T, typename I>
void permute_csx_dny<T,I>::set_matrix_dst(MatrixDenseView_new<T> * M_dst_)
{
    if(M_dst != nullptr) eslog::error("matrix M_dst is already set\n");

    M_dst = M_dst_;
}



template<typename T, typename I>
void permute_csx_dny<T,I>::set_perm_rows(PermutationView_new<I> * perm_rows_)
{
    if(perm_rows != nullptr) eslog::error("perm_rows is already set\n");

    perm_rows = perm_rows_;
}



template<typename T, typename I>
void permute_csx_dny<T,I>::set_perm_cols(PermutationView_new<I> * perm_cols_)
{
    if(perm_cols != nullptr) eslog::error("perm_cols is already set\n");

    perm_cols = perm_cols_;
}



template<typename T, typename I>
void permute_csx_dny<T,I>::perform_zerofill()
{
    stacktimer::push("permute_csx_dny::perform_zerofill");

    if(M_src == nullptr) eslog::error("source matrix is not set\n");
    if(M_dst == nullptr) eslog::error("destination matrix is not set\n");
    if(!M_src->ator->is_data_accessible_cpu()) eslog::error("source matrix must be cpu-accessible\n");
    if(!M_dst->ator->is_data_accessible_cpu()) eslog::error("destination matrix must be cpu-accessible\n");
    if(M_src->nrows != M_dst->nrows || M_src->ncols != M_dst->ncols) eslog::error("matrix sizes dont match\n");
    if(M_src->prop.uplo != M_dst->prop.uplo) eslog::error("matrix uplo does not match\n");
    if(perm_rows != nullptr && perm_rows->size != M_src->nrows) eslog::error("wrong row perm size\n");
    if(perm_cols != nullptr && perm_cols->size != M_src->ncols) eslog::error("wrong col perm size\n");

    fill_dnx<T>::do_all(M_dst, T{0});

    stacktimer::pop();
}



template<typename T, typename I>
void permute_csx_dny<T,I>::perform_copyvals()
{
    stacktimer::push("permute_csx_dny::perform_copyvals");

    if(M_src == nullptr) eslog::error("source matrix is not set\n");
    if(M_dst == nullptr) eslog::error("destination matrix is not set\n");
    if(!M_src->ator->is_data_accessible_cpu()) eslog::error("source matrix must be cpu-accessible\n");
    if(!M_dst->ator->is_data_accessible_cpu()) eslog::error("destination matrix must be cpu-accessible\n");
    if(M_src->nrows != M_dst->nrows || M_src->ncols != M_dst->ncols) eslog::error("matrix sizes dont match\n");
    if(M_src->prop.uplo != M_dst->prop.uplo) eslog::error("matrix uplo does not match\n");
    if(perm_rows != nullptr && perm_rows->size != M_src->nrows) eslog::error("wrong row perm size\n");
    if(perm_cols != nullptr && perm_cols->size != M_src->ncols) eslog::error("wrong col perm size\n");

    bool is_uplo = (M_src->prop.uplo == 'L' || M_src->prop.uplo == 'U');
    bool toggle = ((M_src->order == 'R') == (M_dst->prop.uplo == 'L'));
    size_t dst_ld = M_dst->ld;

    PermutationView_new<I> * perm_primary = (M_src->order == 'R') ? perm_rows : perm_cols;
    PermutationView_new<I> * perm_secdary = (M_src->order == 'R') ? perm_cols : perm_rows;

    eslog::error("not tested, be careful\n");

    for(size_t ips = 0; ips < M_src->get_size_primary(); ips++) {
        I ipd = (perm_primary == nullptr) ? ips : perm_primary->src_to_dst[ips];
        I start = M_src->ptrs[ips];
        I end = M_src->ptrs[ips+1];
        for(I i = start; i < end; i++) {
            I iss = M_src->idxs[i];
            I isd = (perm_secdary == nullptr) ? iss : perm_secdary->src_to_dst[iss];
            T v = M_src->vals[i];

            size_t dst_idx;
            if(!is_uplo || (isd < ipd) == toggle) dst_idx = ipd * dst_ld + isd;
            else dst_idx = isd * dst_ld + ipd;
            M_dst->vals[dst_idx] = v;
        }
    }

    stacktimer::pop();
}



template<typename T, typename I>
void permute_csx_dny<T,I>::perform_all()
{
    perform_zerofill();
    perform_copyvals();
}



template<typename T, typename I>
void permute_csx_dny<T,I>::do_all(MatrixCsxView_new<T,I> * M_src, MatrixDenseView_new<T> * M_dst, PermutationView_new<I> * perm_rows, PermutationView_new<I> * perm_cols)
{
    permute_csx_dny<T,I> instance;
    instance.set_matrix_src(M_src);
    instance.set_matrix_dst(M_dst);
    instance.set_perm_rows(perm_rows);
    instance.set_perm_cols(perm_cols);
    instance.perform_all();
}



#define INSTANTIATE_T_I(T,I) \
template class permute_csx_dny<T,I>;

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
