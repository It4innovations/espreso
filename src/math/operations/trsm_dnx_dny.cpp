
#include "math/operations/trsm_dnx_dny.h"

#include "math/wrappers/math.blas.h"
#include "basis/utilities/stacktimer.h"



namespace espreso {
namespace math {
namespace operations {



template<typename T>
void trsm_dnx_dny<T>::set_system_matrix(MatrixDenseView_new<T> * A_)
{
    if(A != nullptr) eslog::error("matrix A is already set\n");

    A = A_;
}



template<typename T>
void trsm_dnx_dny<T>::set_rhs_sol(MatrixDenseView_new<T> * X_)
{
    if(X != nullptr) eslog::error("matrix X is already set\n");

    X = X_;
}



template<typename T>
void trsm_dnx_dny<T>::perform()
{
    stacktimer::push("trsm_dnx_dny::perform");

    if(A == nullptr) eslog::error("system matrix is not set\n");
    if(X == nullptr) eslog::error("rhs/sol matrix is not set\n");
    if(!A->ator->is_data_accessible_cpu()) eslog::error("matrix A must be cpu-accessible\n");
    if(!X->ator->is_data_accessible_cpu()) eslog::error("matrix X must be cpu-accessible\n");
    if(A->nrows != A->ncols) eslog::error("system matrix has to be square\n");
    if(X->nrows != A->nrows) eslog::error("matrices are incompatible\n");
    if(A->prop.uplo != 'U' && A->prop.uplo != 'L') eslog::error("invalid A uplo\n");
    if(A->prop.diag != 'U' && A->prop.diag != 'N') eslog::error("invalid A diag\n");

    math::blas::trsm(*A, *X);

    stacktimer::pop();
}



template<typename T>
void trsm_dnx_dny<T>::do_all(MatrixDenseView_new<T> * A, MatrixDenseView_new<T> * X)
{
    trsm_dnx_dny<T> instance;
    instance.set_system_matrix(A);
    instance.set_rhs_sol(X);
    instance.perform();
}



#define INSTANTIATE_T(T) \
template class trsm_dnx_dny<T>;

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
