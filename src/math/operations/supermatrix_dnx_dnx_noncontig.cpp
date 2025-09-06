
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
    if(M_src != nullptr) eslog::error("matrix M_src is already set\n");

    M_src = M_src_;
}



template<typename T, typename I>
void supermatrix_dnx_dnx_noncontig<T,I>::set_matrix_destination(MatrixDenseView_new<T> * M_dst_)
{
    if(M_dst != nullptr) eslog::error("matrix M_dst is already set\n");

    M_dst = M_dst_;
}



template<typename T, typename I>
void supermatrix_dnx_dnx_noncontig<T,I>::set_row_map(VectorDenseView_new<I> * row_map_)
{
    if(row_map != nullptr) eslog::error("row_map is already set\n");

    row_map = row_map_;
}



template<typename T, typename I>
void supermatrix_dnx_dnx_noncontig<T,I>::set_col_map(VectorDenseView_new<I> * col_map_)
{
    if(col_map != nullptr) eslog::error("col_map is already set\n");

    col_map = col_map_;
}



template<typename T, typename I>
void supermatrix_dnx_dnx_noncontig<T,I>::set_mode(mode mode_val_)
{
    mode_val = mode_val_;
}



template<typename T, typename I, typename OPERATION>
static void perform_worker(size_t size_primary, size_t size_secdary, const T * __restrict__ src_vals, T * __restrict__ dst_vals, size_t src_ld, size_t dst_ld, const I * primary_map, const I * secdary_map, const OPERATION & op)
{
    if(primary_map == nullptr && secdary_map == nullptr) {
        for(size_t ip = 0; ip < size_primary; ip++) {
            const T * src_ptr = src_vals + ip * src_ld;
            T * dst_ptr = dst_vals + ip * dst_ld;
            for(size_t is = 0; is < size_secdary; is++) {
                op(dst_ptr[is], src_ptr[is]);
            }
        }
    }
    if(primary_map != nullptr && secdary_map == nullptr) {
        for(size_t ips = 0; ips < size_primary; ips++) {
            size_t ipd = primary_map[ips];
            const T * src_ptr = src_vals + ips * src_ld;
            T * dst_ptr = dst_vals + ipd * dst_ld;
            for(size_t is = 0; is < size_secdary; is++) {
                op(dst_ptr[is], src_ptr[is]);
            }
        }
    }
    if(primary_map == nullptr && secdary_map != nullptr) {
        for(size_t ip = 0; ip < size_primary; ip++) {
            const T * src_ptr = src_vals + ip * src_ld;
            T * dst_ptr = dst_vals + ip * dst_ld;
            for(size_t iss = 0; iss < size_secdary; iss++) {
                size_t isd = secdary_map[iss];
                op(dst_ptr[isd], src_ptr[iss]);
            }
        }
    }
    if(primary_map != nullptr && secdary_map != nullptr) {
        for(size_t ips = 0; ips < size_primary; ips++) {
            size_t ipd = primary_map[ips];
            const T * src_ptr = src_vals + ips * src_ld;
            T * dst_ptr = dst_vals + ipd * dst_ld;
            for(size_t iss = 0; iss < size_secdary; iss++) {
                size_t isd = secdary_map[iss];
                op(dst_ptr[isd], src_ptr[iss]);
            }
        }
    }
}



template<typename T, typename I>
void supermatrix_dnx_dnx_noncontig<T,I>::perform()
{
    stacktimer::push("supermatrix_dnx_dnx_noncontig::perform");

    if(M_src == nullptr) eslog::error("source matrix is not set\n");
    if(M_dst == nullptr) eslog::error("destination matrix is not set\n");
    if(!M_src->ator->is_data_accessible_cpu()) eslog::error("source matrix must be cpu-accessible\n");
    if(!M_dst->ator->is_data_accessible_cpu()) eslog::error("destination matrix must be cpu-accessible\n");
    if(row_map != nullptr && !row_map->ator->is_data_accessible_cpu()) eslog::error("row_map must be cpu-accessible\n");
    if(col_map != nullptr && !col_map->ator->is_data_accessible_cpu()) eslog::error("col_map must be cpu-accessible\n");
    if(M_src->order != M_dst->order) eslog::error("matrix orders do not match\n");
    size_t nrows = ((row_map == nullptr) ? M_dst->nrows : row_map->size);
    size_t ncols = ((col_map == nullptr) ? M_dst->ncols : col_map->size);
    if(M_src->nrows != nrows || M_src->ncols != ncols) eslog::error("wrong source matrix size\n");

    VectorDenseView_new<I> * primary_map = ((M_dst->order == 'R') ? row_map : col_map);
    VectorDenseView_new<I> * secdary_map = ((M_dst->order == 'R') ? col_map : row_map);
    I * primary_map_vals = ((primary_map != nullptr) ? primary_map->vals : nullptr);
    I * secdary_map_vals = ((secdary_map != nullptr) ? secdary_map->vals : nullptr);

    switch(mode_val) {
        case mode::assign:
            perform_worker(M_src->get_size_primary(), M_src->get_size_secdary(), M_src->vals, M_dst->vals, M_src->ld, M_dst->ld, primary_map_vals, secdary_map_vals, [](T & dst, T src) { dst = src; });
            break;
        case mode::accumulate:
            perform_worker(M_src->get_size_primary(), M_src->get_size_secdary(), M_src->vals, M_dst->vals, M_src->ld, M_dst->ld, primary_map_vals, secdary_map_vals, [](T & dst, T src) { dst += src; });
            break;
        case mode::accumulate_atomic:
            perform_worker(M_src->get_size_primary(), M_src->get_size_secdary(), M_src->vals, M_dst->vals, M_src->ld, M_dst->ld, primary_map_vals, secdary_map_vals, [](T & dst, T src) { utils::atomic_add(dst, src); });
            break;
        default:
            eslog::error("wrong mode\n");
    }

    stacktimer::pop();
}



template<typename T, typename I>
void supermatrix_dnx_dnx_noncontig<T,I>::do_all(MatrixDenseView_new<T> * M_src, MatrixDenseView_new<T> * M_dst, VectorDenseView_new<I> * row_map, VectorDenseView_new<I> * col_map, mode mode_val)
{
    supermatrix_dnx_dnx_noncontig<T,I> instance;
    instance.set_matrix_source(M_src);
    instance.set_matrix_destination(M_dst);
    instance.set_row_map(row_map);
    instance.set_col_map(col_map);
    instance.set_mode(mode_val);
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
