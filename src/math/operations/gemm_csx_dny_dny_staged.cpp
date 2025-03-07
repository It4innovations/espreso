
#include "math/operations/gemm_csx_dny_dny_staged.h"

#include "basis/utilities/stacktimer.h"



namespace espreso {
namespace math {
namespace operations {



template<typename T, typename I>
gemm_csx_dny_dny_staged<T,I>::~gemm_csx_dny_dny_staged()
{
    finalize();
}



template<typename T, typename I>
void gemm_csx_dny_dny_staged<T,I>::set_matrix_A(MatrixCsxView_new<T,I> * A_)
{
    A = A_;
}



template<typename T, typename I>
void gemm_csx_dny_dny_staged<T,I>::set_matrix_B(MatrixDenseView_new<T> * B_)
{
    B = B_;
}



template<typename T, typename I>
void gemm_csx_dny_dny_staged<T,I>::set_matrix_C(MatrixDenseView_new<T> * C_)
{
    C = C_;
}



template<typename T, typename I>
void gemm_csx_dny_dny_staged<T,I>::set_coefficients(T alpha_, T beta_)
{
    alpha = alpha_;
    beta = beta_;
}



template<typename T, typename I>
void gemm_csx_dny_dny_staged<T,I>::preprocess()
{
    if(preprocess_called) eslog::error("preproces has already been called\n");
    if(A == nullptr) eslog::error("matrix A is not set\n");
    if(B == nullptr) eslog::error("matrix B is not set\n");
    if(C == nullptr) eslog::error("matrix C is not set\n");
    if(A->nrows != C->nrows || B->ncols != C->ncols || A->ncols != B->nrows) eslog::error("incompatible matrix sizes\n");
    if(B->order != C->order) eslog::error("order of B and C must match\n");

    stacktimer::push("gemm_csx_dny_dny_staged::preprocess");

    spblas::mm(*A, *B, *C, alpha, beta, handle_abc, 'P');

    stacktimer::pop();

    preprocess_called = true;
}



template<typename T, typename I>
void gemm_csx_dny_dny_staged<T,I>::perform()
{
    if(!preprocess_called) eslog::error("preprocess has not been called\n");
    if(A == nullptr) eslog::error("matrix A is not set\n");
    if(B == nullptr) eslog::error("matrix B is not set\n");
    if(C == nullptr) eslog::error("matrix C is not set\n");
    if(A->nrows != C->nrows || B->ncols != C->ncols || A->ncols != B->nrows) eslog::error("incompatible matrix sizes\n");
    if(B->order != C->order) eslog::error("order of B and C must match\n");

    stacktimer::push("gemm_csx_dny_dny_staged::perform");

    spblas::mm(*A, *B, *C, alpha, beta, handle_abc, 'C');

    stacktimer::pop();
}



template<typename T, typename I>
void gemm_csx_dny_dny_staged<T,I>::finalize()
{
    if(preprocess_called) {
        spblas::mm(*A, *B, *C, alpha, beta, handle_abc, 'F');
    }
    preprocess_called = false;
}



#define INSTANTIATE_T_I(T,I) \
template class gemm_csx_dny_dny_staged<T,I>;

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
