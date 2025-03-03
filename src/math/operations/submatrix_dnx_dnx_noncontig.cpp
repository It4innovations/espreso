
#include "math/operations/submatrix_dnx_dnx_noncontig.h"

#include <algorithm>

#include "math/operations/copy_dnx.h"



namespace espreso {
namespace math {
namespace operations {



template<typename T, typename I>
void submatrix_dnx_dnx_noncontig<T,I>::set_matrix_source(MatrixDenseView_new<T> * M_src_)
{
    M_src = M_src_;
}



template<typename T, typename I>
void submatrix_dnx_dnx_noncontig<T,I>::set_matrix_destination(MatrixDenseView_new<T> * M_dst_)
{
    M_dst = M_dst_;
}



template<typename T, typename I>
void submatrix_dnx_dnx_noncontig<T,I>::set_row_map(VectorDenseView_new<I> * row_map_)
{
    row_map = row_map_;
}



template<typename T, typename I>
void submatrix_dnx_dnx_noncontig<T,I>::set_col_map(VectorDenseView_new<I> * col_map_)
{
    col_map = col_map_;
}



template<typename T, typename I>
void submatrix_dnx_dnx_noncontig<T,I>::perform()
{
    if(M_src == nullptr) eslog::error("source matrix is not set\n");
    if(M_dst == nullptr) eslog::error("destination matrix is not set\n");
    if(M_src->order != M_dst->order) eslog::error("matrix orders do not match\n");
    size_t nrows = ((row_map == nullptr) ? M_src->nrows : row_map->size);
    size_t ncols = ((col_map == nullptr) ? M_src->ncols : col_map->size);
    if(M_dst->nrows != nrows || M_dst->ncols != ncols) eslog::error("wrong destination matrix size\n");

    size_t dst_size_primary = M_dst->get_size_primary();
    size_t dst_size_secdary = M_dst->get_size_secdary();

    T * src_vals = M_src->vals;
    T * dst_vals = M_dst->vals;
    size_t src_ld = M_src->ld;
    size_t dst_ld = M_dst->ld;

    VectorDenseView_new<I> * primary_map = ((M_src->order == 'R') ? row_map : col_map);
    VectorDenseView_new<I> * secdary_map = ((M_src->order == 'R') ? col_map : row_map);

    if(primary_map == nullptr && secdary_map == nullptr) {
        copy_dnx<T>::do_all(M_src, M_dst);
    }
    if(primary_map != nullptr && secdary_map == nullptr) {
        I * subset_primary = primary_map->vals;
        for(size_t ipd = 0; ipd < dst_size_primary; ipd++) {
            size_t ips = subset_primary[ipd];
            std::copy_n(src_vals + ips * src_ld, dst_size_secdary, dst_vals + ipd * dst_ld);
        }
    }
    if(primary_map == nullptr && secdary_map != nullptr) {
        I * subset_secdary = secdary_map->vals;
        for(size_t ip = 0; ip < dst_size_primary; ip++) {
            for(size_t isd = 0; isd < dst_size_secdary; isd++) {
                size_t iss = subset_secdary[isd];
                dst_vals[ip * dst_ld + isd] = src_vals[ip * src_ld + iss];
            }
        }
    }
    if(primary_map != nullptr && secdary_map != nullptr) {
        I * subset_primary = primary_map->vals;
        I * subset_secdary = secdary_map->vals;
        for(size_t ipd = 0; ipd < dst_size_primary; ipd++) {
            size_t ips = subset_primary[ipd];
            for(size_t isd = 0; isd < dst_size_secdary; isd++) {
                size_t iss = subset_secdary[isd];
                dst_vals[ipd * dst_ld + isd] = src_vals[ips * src_ld + iss];
            }
        }
    }
}



template<typename T, typename I>
void submatrix_dnx_dnx_noncontig<T,I>::do_all(MatrixDenseView_new<T> * M_src, MatrixDenseView_new<T> * M_dst, VectorDenseView_new<I> * row_map, VectorDenseView_new<I> * col_map)
{
    submatrix_dnx_dnx_noncontig<T,I> instance;
    instance.set_matrix_source(M_src);
    instance.set_matrix_destination(M_dst);
    instance.set_row_map(row_map);
    instance.set_col_map(col_map);
    instance.perform();
}



#define INSTANTIATE_T_I(T,I) \
template class submatrix_dnx_dnx_noncontig<T,I>;

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
