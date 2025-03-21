
#include "math/operations/supermatrix_dnx_dnx_noncontig.h"

#include <algorithm>

#include "math/operations/copy_dnx.h"
#include "basis/utilities/stacktimer.h"



namespace espreso {
namespace math {
namespace operations {



template<typename T, typename I>
void supermatrix_dnx_dnx_noncontig<T,I>::set_matrix_source(MatrixDenseView_new<T> * M_src_)
{
    M_src = M_src_;
}



template<typename T, typename I>
void supermatrix_dnx_dnx_noncontig<T,I>::set_matrix_destination(MatrixDenseView_new<T> * M_dst_)
{
    M_dst = M_dst_;
}



template<typename T, typename I>
void supermatrix_dnx_dnx_noncontig<T,I>::set_row_map(VectorDenseView_new<I> * row_map_)
{
    row_map = row_map_;
}



template<typename T, typename I>
void supermatrix_dnx_dnx_noncontig<T,I>::set_col_map(VectorDenseView_new<I> * col_map_)
{
    col_map = col_map_;
}



template<typename T, typename I>
void supermatrix_dnx_dnx_noncontig<T,I>::perform()
{
    if(M_src == nullptr) eslog::error("source matrix is not set\n");
    if(M_dst == nullptr) eslog::error("destination matrix is not set\n");
    if(M_src->order != M_dst->order) eslog::error("matrix orders do not match\n");
    size_t nrows = ((row_map == nullptr) ? M_dst->nrows : row_map->size);
    size_t ncols = ((col_map == nullptr) ? M_dst->ncols : col_map->size);
    if(M_src->nrows != nrows || M_src->ncols != ncols) eslog::error("wrong source matrix size\n");

    stacktimer::push("supermatrix_dnx_dnx_noncontig::perform");

    size_t src_size_primary = M_src->get_size_primary();
    size_t src_size_secdary = M_src->get_size_secdary();

    T * src_vals = M_src->vals;
    T * dst_vals = M_dst->vals;
    size_t src_ld = M_src->ld;
    size_t dst_ld = M_dst->ld;

    VectorDenseView_new<I> * primary_map = ((M_dst->order == 'R') ? row_map : col_map);
    VectorDenseView_new<I> * secdary_map = ((M_dst->order == 'R') ? col_map : row_map);

    if(primary_map == nullptr && secdary_map == nullptr) {
        copy_dnx<T>::do_all(M_src, M_dst);
    }
    if(primary_map != nullptr && secdary_map == nullptr) {
        I * subset_primary = primary_map->vals;
        for(size_t ips = 0; ips < src_size_primary; ips++) {
            size_t ipd = subset_primary[ips];
            std::copy_n(src_vals + ips * src_ld, src_size_secdary, dst_vals + ipd * dst_ld);
        }
    }
    if(primary_map == nullptr && secdary_map != nullptr) {
        I * subset_secdary = secdary_map->vals;
        for(size_t ip = 0; ip < src_size_primary; ip++) {
            for(size_t iss = 0; iss < src_size_secdary; iss++) {
                size_t isd = subset_secdary[iss];
                dst_vals[ip * dst_ld + isd] = src_vals[ip * src_ld + iss];
            }
        }
    }
    if(primary_map != nullptr && secdary_map != nullptr) {
        I * subset_primary = primary_map->vals;
        I * subset_secdary = secdary_map->vals;
        for(size_t ips = 0; ips < src_size_primary; ips++) {
            size_t ipd = subset_primary[ips];
            for(size_t iss = 0; iss < src_size_secdary; iss++) {
                size_t isd = subset_secdary[iss];
                dst_vals[ipd * dst_ld + isd] = src_vals[ips * src_ld + iss];
            }
        }
    }

    stacktimer::pop();
}



template<typename T, typename I>
void supermatrix_dnx_dnx_noncontig<T,I>::do_all(MatrixDenseView_new<T> * M_src, MatrixDenseView_new<T> * M_dst, VectorDenseView_new<I> * row_map, VectorDenseView_new<I> * col_map)
{
    supermatrix_dnx_dnx_noncontig<T,I> instance;
    instance.set_matrix_source(M_src);
    instance.set_matrix_destination(M_dst);
    instance.set_row_map(row_map);
    instance.set_col_map(col_map);
    instance.perform();
}



#define INSTANTIATE_T_I(T,I) \
template class supermatrix_dnx_dnx_noncontig<T,I>;

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
