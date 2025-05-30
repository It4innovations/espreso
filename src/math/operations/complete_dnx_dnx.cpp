
#include "math/operations/complete_dnx_dnx.h"

#include "math/primitives_new/allocator_new.h"
#include "math/primitives_new/matrix_dense_data_new.h"
#include "math/operations/convert_dnx_dny.h"
#include "math/operations/copy_dnx.h"
#include "math/operations/transpose_dnx_dnx.h"
#include "basis/utilities/stacktimer.h"



namespace espreso {
namespace math {
namespace operations {



template<typename T>
void complete_dnx_dnx<T>::set_matrix_src(MatrixDenseView_new<T> * M_src_)
{
    if(M_src != nullptr) eslog::error("source matrix is already set\n");

    M_src = M_src_;
}



template<typename T>
void complete_dnx_dnx<T>::set_matrix_dst(MatrixDenseView_new<T> * M_dst_)
{
    if(M_dst != nullptr) eslog::error("destination matrix is already set\n");

    M_dst = M_dst_;
}



template<typename T>
void complete_dnx_dnx<T>::set_conj(bool do_conj_)
{
    do_conj = do_conj_;
}



template<typename T>
void complete_dnx_dnx<T>::perform()
{
    stacktimer::push("complete_dnx_dnx::perform");

    if(M_src == nullptr) eslog::error("source matrix is not set\n");
    if(M_dst == nullptr) eslog::error("destination matrix is not set\n");
    if(M_src == M_dst) eslog::error("in-place is not supported\n");
    if(!M_src->ator->is_data_accessible_cpu()) eslog::error("source matrix must be cpu-accessible\n");
    if(!M_dst->ator->is_data_accessible_cpu()) eslog::error("destination matrix must be cpu-accessible\n");
    if(M_src->nrows != M_dst->nrows || M_src->ncols != M_dst->ncols) eslog::error("matrix sizes dont match\n");
    if(M_src->nrows != M_src->ncols) eslog::error("matrices must be square\n");
    if(M_src->order != M_dst->order) eslog::error("matrix order does not match\n");
    if(M_src->prop.uplo != 'L' && M_src->prop.uplo != 'U') eslog::error("source matrix must be triangular\n");

    transpose_dnx_dnx<T>::do_all(M_src, M_dst, do_conj);

    MatrixDenseView_new<T> M_dst_2 = *M_dst;
    M_dst_2.prop.uplo = M_src->prop.uplo;

    copy_dnx<T>::do_all(M_src, &M_dst_2, false);

    stacktimer::pop();
}



template<typename T>
void complete_dnx_dnx<T>::do_all(MatrixDenseView_new<T> * M_src, MatrixDenseView_new<T> * M_dst, bool do_conj)
{
    complete_dnx_dnx<T> instance;
    instance.set_matrix_src(M_src);
    instance.set_matrix_dst(M_dst);
    instance.set_conj(do_conj);
    instance.perform();
}



#define INSTANTIATE_T(T) \
template class complete_dnx_dnx<T>;

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

