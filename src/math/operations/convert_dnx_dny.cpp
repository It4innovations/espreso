
#include "math/operations/convert_dnx_dny.h"

#include "math/operations/copy_dnx.h"
#include "math/wrappers/math.blas.h"
#include "basis/utilities/stacktimer.h"



namespace espreso {
namespace math {
namespace operations {



template<typename T>
void convert_dnx_dny<T>::set_matrix_src(MatrixDenseView_new<T> * M_src_)
{
    M_src = M_src_;
}



template<typename T>
void convert_dnx_dny<T>::set_matrix_dst(MatrixDenseView_new<T> * M_dst_)
{
    M_dst = M_dst_;
}



template<typename T>
void convert_dnx_dny<T>::set_conj(bool do_conj_)
{
    do_conj = do_conj_;
}



template<typename T>
void convert_dnx_dny<T>::perform()
{
    if(M_src == nullptr) eslog::error("source matrix is not set\n");
    if(M_dst == nullptr) eslog::error("destination matrix is not set\n");
    if(M_src->nrows != M_dst->nrows || M_src->ncols != M_dst->ncols) eslog::error("matrix sizes dont match\n");

    stacktimer::push("convert_dnx_dny::perform");

    if(M_src->order == M_dst->order) {
        copy_dnx<T>::do_all(M_src, M_dst, do_conj);
    }
    else {
        blas::transpose(M_src->nrows, M_src->ncols, M_src->vals, M_src->ld, M_dst->vals, M_dst->ld, M_src->order, do_conj);
    }

    stacktimer::pop();
}



template<typename T>
void convert_dnx_dny<T>::do_all(MatrixDenseView_new<T> * M_src, MatrixDenseView_new<T> * M_dst, bool do_conj)
{
    convert_dnx_dny<T> instance;
    instance.set_matrix_src(M_src);
    instance.set_matrix_dst(M_dst);
    instance.set_conj(do_conj);
    instance.perform();
}



#define INSTANTIATE_T(T) \
template class convert_dnx_dny<T>;

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
