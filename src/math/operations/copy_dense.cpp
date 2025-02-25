
#include "math/operations/copy_dense.h"



namespace espreso {
namespace math {
namespace operations {



template<typename T>
void copy_dense<T>::set_matrix_src(MatrixDenseView_new<T> * M_src_)
{
    M_src = M_src_;
}



template<typename T>
void copy_dense<T>::set_matrix_dst(MatrixDenseView_new<T> * M_dst_)
{
    M_dst = M_dst_;
}



template<typename T>
void copy_dense<T>::perform()
{
    if(M_src->nrows != M_dst->nrows || M_src->ncols != M_dst->ncols) eslog::error("matrix sizes dont match\n");
    if(M_src->order != M_dst->order) eslog::error("matrix order dont match\n");

    size_t num_blocks = M_src->get_num_blocks();
    size_t block_size = M_src->get_block_size();
    for(size_t i = 0; i < num_blocks; i++) {
        std::copy_n(M_src->vals + i * M_src->ld, block_size, M_dst->vals + i * M_dst->ld);
    }
}



template<typename T>
void copy_dense<T>::do_all(MatrixDenseView_new<T> * M_src, MatrixDenseView_new<T> * M_dst)
{
    copy_dense<T> instance;
    instance.set_src(M_src);
    instance.set_dst(M_dst);
    instance.perform();
}



#define INSTANTIATE_T(T) \
template class copy_dense<T>;

    #define INSTANTIATE \
    /* INSTANTIATE_T(float) */ \
    INSTANTIATE_T(double) \
    /* INSTANTIATE_T(std::complex<float>) */ \
    /* INSTANTIATE_T(std::complex<double>) */

        INSTANTIATE

    #undef INSTANTIATE
#undef INSTANTIATE_T



}
}
}
