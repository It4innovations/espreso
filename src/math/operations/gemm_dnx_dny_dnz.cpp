
#include "math/operations/gemm_dnx_dny_dnz.h"

#include "math/wrappers/math.blas.h"
#include "basis/utilities/stacktimer.h"



namespace espreso {
namespace math {
namespace operations {



template<typename T>
void gemm_dnx_dny_dnz<T>::set_matrix_A(MatrixDenseView_new<T> * A_)
{
    if(A != nullptr) eslog::error("matrix A is already set\n");

    A = A_;
}



template<typename T>
void gemm_dnx_dny_dnz<T>::set_matrix_B(MatrixDenseView_new<T> * B_)
{
    if(B != nullptr) eslog::error("matrix B is already set\n");

    B = B_;
}



template<typename T>
void gemm_dnx_dny_dnz<T>::set_matrix_C(MatrixDenseView_new<T> * C_)
{
    if(C != nullptr) eslog::error("matrix C is already set\n");

    C = C_;
}



template<typename T>
void gemm_dnx_dny_dnz<T>::set_coefficients(T alpha_, T beta_)
{
    alpha = alpha_;
    beta = beta_;
}



template<typename T>
void gemm_dnx_dny_dnz<T>::set_conj(bool conj_A_, bool conj_B_)
{
    conj_A = conj_A_;
    conj_B = conj_B_;
}



template<typename T>
void gemm_dnx_dny_dnz<T>::perform()
{
    stacktimer::push("gemm_dnx_dny_dnz::perform");

    if(A == nullptr) eslog::error("matrix A is not set\n");
    if(B == nullptr) eslog::error("matrix B is not set\n");
    if(C == nullptr) eslog::error("matrix C is not set\n");
    if(!A->ator->is_data_accessible_cpu()) eslog::error("matrix A must be cpu-accessible\n");
    if(!B->ator->is_data_accessible_cpu()) eslog::error("matrix B must be cpu-accessible\n");
    if(!C->ator->is_data_accessible_cpu()) eslog::error("matrix C must be cpu-accessible\n");
    if(A->nrows != C->nrows || B->ncols != C->ncols || A->ncols != B->nrows) eslog::error("incompatible matrices");
    if(utils::is_complex<T>() && (conj_A || conj_B)) eslog::error("conjugation not supported yet\n");

    math::blas::gemm(*A, *B, *C, alpha, beta);

    stacktimer::pop();
}



template<typename T>
void gemm_dnx_dny_dnz<T>::do_all(MatrixDenseView_new<T> * A, MatrixDenseView_new<T> * B, MatrixDenseView_new<T> * C, T alpha, T beta, bool conj_A, bool conj_B)
{
    gemm_dnx_dny_dnz<T> instance;
    instance.set_matrix_A(A);
    instance.set_matrix_B(B);
    instance.set_matrix_C(C);
    instance.set_coefficients(alpha, beta);
    instance.set_conj(conj_A, conj_B);
    instance.perform();
}



#define INSTANTIATE_T(T) \
template class gemm_dnx_dny_dnz<T>;

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
