
#include "math/operations/complete_dnx.h"



namespace espreso {
namespace math {
namespace operations {



template<typename T>
void complete_dnx<T>::set_matrix(MatrixDenseView_new<T> * M_)
{
    M = M_;
}



template<typename T>
void complete_dnx<T>::set_conj(bool do_conj_)
{
    do_conj = do_conj_;
}



template<typename T>
void complete_dnx<T>::perform()
{
    if(M == nullptr) eslog::error("matrix is not set\n");
    if(M->nrows != M->ncols) eslog::error("matrix must be square\n");
    if(M->uplo != 'L' && M->uplo != 'U') eslog::error("matrix must be upper or lower\n");

    size_t size = M->nrows;

    MatrixDenseData_new<T> M2;
    M2.set(size, size, change_order(M->order), AllocatorCPU_new::get_singleton());
    M2.alloc();

    convert_dnx_dny<T>::do_all(M, &M2, do_conj);

    M2.prop.uplo = change_uplo(M->uplo);

    copy_dnx<T>::do_all(&M2, M, false);

    M2.clear();
}



template<typename T>
void complete_dnx<T>::do_all(MatrixDenseView_new<T> * M, bool do_conj)
{
    complete_dnx<T> instance;
    instance.set_matrix(M);
    instance.set_conj(do_conj);
    instance.perform();
}



#define INSTANTIATE_T(T) \
template class complete_dnx<T>;

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

