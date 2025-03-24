
#include "math/operations/herk_dnx_dny.h"

#include "basis/utilities/stacktimer.h"



namespace espreso {
namespace math {
namespace operations {



template<typename T>
void herk_dnx_dny<T>::set_matrix_A(MatrixDenseView_new<T> * A_)
{
    A = A_;
}



template<typename T>
void herk_dnx_dny<T>::set_matrix_C(MatrixDenseView_new<T> * C_)
{
    C = C_;
}



template<typename T>
void herk_dnx_dny<T>::set_mode(blas::herk_mode mode_)
{
    mode = mode_;
    mode_set = true;
}



template<typename T>
void herk_dnx_dny<T>::set_coefficients(Treal alpha_, Treal beta_)
{
    alpha = alpha_;
    beta = beta_;
}



template<typename T>
void herk_dnx_dny<T>::perform()
{
    stacktimer::push("herk_dnx_dny::perform");

    if(A == nullptr) eslog::error("matrix A is not set\n");
    if(C == nullptr) eslog::error("matrix C is not set\n");
    if(!mode_set) eslog::error("mode is not set");
    if(C->nrows != C->ncols) eslog::error("C must be square\n");
    if(mode == blas::herk_mode::AAh && A->nrows != C->nrows) eslog::error("incompatible matrix sizes\n");
    if(mode == blas::herk_mode::AhA && A->ncols != C->ncols) eslog::error("incompatible matrix sizes\n");
    if(C->prop.uplo != 'U' && C->prop.uplo != 'L') eslog::error("C uplo is not set\n");

    blas::herk(*A, *C, mode, alpha, beta);

    stacktimer::pop();
}



template<typename T>
void herk_dnx_dny<T>::do_all(MatrixDenseView_new<T> * A, MatrixDenseView_new<T> * C, blas::herk_mode mode, Treal alpha, Treal beta)
{
    herk_dnx_dny<T> instance;
    instance.set_matrix_A(A);
    instance.set_matrix_C(C);
    instance.set_mode(mode);
    instance.set_coefficients(alpha, beta);
    instance.perform();
}



#define INSTANTIATE_T(T) \
template class herk_dnx_dny<T>;

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
