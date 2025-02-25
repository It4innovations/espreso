
#include "math/operations/copy_csx.h"



namespace espreso {
namespace math {
namespace operations {



template<typename T, typename I>
void copy_csx<T,I>::set_matrix_src(MatrixCsxView_new<T,I> * M_src_)
{
    M_src = M_src_;
}



template<typename T, typename I>
void copy_csx<T,I>::set_matrix_dst(MatrixCsxView_new<T,I> * M_dst_)
{
    M_dst = M_dst;
}



template<typename T, typename I>
void copy_csx<T,I>::perform()
{
    if(M_src == nullptr) eslog::error("source matrix is not set\n");
    if(M_dst == nullptr) eslog::error("destination matrix is not set\n");
    if(M_src->nrows != M_dst.nrows || M_src.ncols != M_dst.ncols || M_src.nnz != M_dst.nnz) eslog::error("matrix sizes dont match\n");
    if(M_src->order != M_dst.order) eslog::error("matrix orders dont match\n");

    std::copy_n(M_src->ptrs, M_src->get_primary_size() + 1, M_dst->ptrs);
    std::copy_n(M_src->idxs, M_src->nnz, M_dst->idxs);
    std::copy_n(M_src->vals, M_src->nnz, M_dst->vals);
}



template<typename T, typename I>
void copy_csx<T,I>::do_all(MatrixCsxView_new<T,I> * M_src, MatrixCsxView_new<T,I> * M_dst)
{
    copy_csx<T,I> instance;
    instance.set_matrix_src(M_src);
    instance.set_matrix_dst(M_dst);
    instance.perform();
}



#define INSTANTIATE_T_I(T,I) \
template class copy_csx<T,I>;

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
