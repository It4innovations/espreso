
#include "math/operations/trsm_csx_dny_staged.h"

#include "math/operations/copy_dnx.h"



namespace espreso {
namespace math {
namespace operations {



template<typename T, typename I>
trsm_csx_dny_staged<T,I>::~trsm_csx_dny_staged()
{
    finalize();
}



template<typename T, typename I>
void trsm_csx_dny_staged<T,I>::set_system_matrix(MatrixCsxView_new<T,I> * A_)
{
    A = A_;
}



template<typename T, typename I>
void trsm_csx_dny_staged<T,I>::set_rhs_sol(MatrixDenseView_new<T> * X_)
{
    X = X_;
}



template<typename T, typename I>
void trsm_csx_dny_staged<T,I>::preprocess()
{
    if(A == nullptr) eslog::error("matrix A is not set\n");
    if(X == nullptr) eslog::error("matrix X is not set\n");
    if(preprocess_called) eslog::error("preprocess has already been called\n");
    if(A->prop.uplo != 'U' && A->prop.uplo != 'L') eslog::error("invalid A uplo\n");
    if(A->prop.diag != 'U' && A->prop.diag != 'N') eslog::error("invalid A diag\n");

    Y.set(X->nrows, X->ncols, X->order, AllocatorCPU_new::get_singleton());

    math::spblas::trsm(*A, *X, Y, handle, 'P');

    preprocess_called = true;
}



template<typename T, typename I>
void trsm_csx_dny_staged<T,I>::perform()
{
    if(A == nullptr) eslog::error("matrix A is not set\n");
    if(X == nullptr) eslog::error("matrix X is not set\n");
    if(!preprocess_called) eslog::error("preprocess has not been called\n");
    if(A->prop.uplo != 'U' && A->prop.uplo != 'L') eslog::error("invalid A uplo\n");
    if(A->prop.diag != 'U' && A->prop.diag != 'N') eslog::error("invalid A diag\n");

    Y.alloc();

    math::spblas::trsm(*A, *X, Y, handle, 'C');

    copy_dnx<T>::do_all(X, &Y);

    Y.free();
}



template<typename T, typename I>
void trsm_csx_dny_staged<T,I>::finalize()
{
    if(preprocess_called) {
        math::spblas::trsm(*A, *X, Y, handle, 'F');
    }
    preprocess_called = false;
}



#define INSTANTIATE_T_I(T,I) \
template class trsm_csx_dny_staged<T,I>;

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
