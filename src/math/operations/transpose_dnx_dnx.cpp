
#include "math/operations/transpose_dnx_dnx.h"
#include "math/wrappers/math.blas.h"
#include "basis/utilities/stacktimer.h"



namespace espreso {
namespace math {
namespace operations {



template<typename T>
void transpose_dnx_dnx<T>::set_matrix_src(MatrixDenseView_new<T> * M_src_)
{
    if(M_src != nullptr) eslog::error("source matrix is already set\n");

    M_src = M_src_;
}



template<typename T>
void transpose_dnx_dnx<T>::set_matrix_dst(MatrixDenseView_new<T> * M_dst_)
{
    if(M_dst != nullptr) eslog::error("destination matrix is already set\n");

    M_dst = M_dst_;
}



template<typename T>
void transpose_dnx_dnx<T>::set_conj(bool do_conj_)
{
    do_conj = do_conj_;
}



template<typename T>
void transpose_dnx_dnx<T>::perform()
{
    stacktimer::push("transpose_dnx_dnx::perform");

    if(M_src == nullptr) eslog::error("source matrix is not set\n");
    if(M_dst == nullptr) eslog::error("destination matrix is not set\n");
    if(!M_src->ator->is_data_accessible_cpu()) eslog::error("source matrix must be cpu-accessible\n");
    if(!M_dst->ator->is_data_accessible_cpu()) eslog::error("destination matrix must be cpu-accessible\n");
    if(M_src->nrows != M_dst->ncols || M_src->ncols != M_dst->nrows) eslog::error("incompatible matrix sizes\n");
    if(M_src->order != M_dst->order) eslog::error("matrices must have equal order\n"); // if different order, transposet reordered view is enough

    if(M_src == M_dst) {
        if(M_src->nrows != M_src->ncols) eslog::error("in-place works only for square matrices\n");

        math::blas::transpose_inplace(M_src->nrows, M_src->vals, M_src->ld, M_src->order, do_conj);
    }
    else {
        math::blas::transpose(M_src->nrows, M_src->ncols, M_src->vals, M_src->ld, M_dst->vals, M_dst->ld, M_src->order, do_conj);
    }

    stacktimer::pop();
}



template<typename T>
void transpose_dnx_dnx<T>::do_all(MatrixDenseView_new<T> * M_src, MatrixDenseView_new<T> * M_dst, bool do_conj)
{
    transpose_dnx_dnx<T> instance;
    instance.set_matrix_src(M_src);
    instance.set_matrix_dst(M_dst);
    instance.set_conj(do_conj);
    instance.perform();
}



#define INSTANTIATE_T(T) \
template class transpose_dnx_dnx<T>;

    #define INSTANTIATE \
    /* INSTANTIATE_T(float) */ \
    INSTANTIATE_T(double) \
    /* INSTANTIATE_T(std::complex<float>) */ \
    INSTANTIATE_T(std::complex<double>)

        INSTANTIATE

    #undef INSTANTIATE
#undef INSTANTIATE_T



}
}
}
