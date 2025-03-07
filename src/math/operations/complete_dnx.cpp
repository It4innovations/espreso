
#include "math/operations/complete_dnx.h"

#include "math/primitives_new/matrix_dense_data_new.h"
#include "math/operations/convert_dnx_dny.h"
#include "math/operations/copy_dnx.h"
#include "basis/utilities/stacktimer.h"



namespace espreso {
namespace math {
namespace operations {



template<typename T>
void complete_dnx<T>::set_matrix(MatrixDenseView_new<T> * M_)
{
    M = M_;
}



template<typename T>
void complete_dnx<T>::set_orig_uplo(char orig_uplo_)
{
    orig_uplo = orig_uplo_;
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
    if(orig_uplo != 'L' && orig_uplo != 'U') eslog::error("original uplo is not set\n");

    stacktimer::push("complete_dnx::perform");

    size_t size = M->nrows;

    MatrixDenseView_new<T> M2 = *M; // same as M, but ignore uplo
    M2.prop.uplo = 'F';

    MatrixDenseData_new<T> M3; // M2 but converted to the other order
    M3.set(size, size, change_order(M2.order), AllocatorCPU_new::get_singleton());
    M3.alloc();
    convert_dnx_dny<T>::do_all(&M2, &M3, do_conj);

    MatrixDenseView_new<T> M4 = M3.get_transposed_reordered_view(); // same order as M, but the uplo is opposite
    M4.prop.uplo = change_uplo(orig_uplo);

    MatrixDenseView_new<T> M5 = *M; // same as M, but has opposite uplo
    M5.prop.uplo = change_uplo(orig_uplo);

    copy_dnx<T>::do_all(&M4, &M5);

    M3.clear();

    stacktimer::pop();
}



template<typename T>
void complete_dnx<T>::do_all(MatrixDenseView_new<T> * M, char orig_uplo, bool do_conj)
{
    complete_dnx<T> instance;
    instance.set_matrix(M);
    instance.set_orig_uplo(orig_uplo);
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

