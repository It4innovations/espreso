
#include "math/operations/gemm_csx_dny_dny_staged.h"

#include "basis/utilities/stacktimer.h"



namespace espreso {
namespace math {
namespace operations {



template<typename T, typename I>
void gemm_csx_dny_dny_staged<T,I>::set_matrix_A(MatrixCsxView_new<T,I> * A_)
{
    if(A != nullptr) eslog::error("matrix A is already set\n");

    A = A_;
}



template<typename T, typename I>
void gemm_csx_dny_dny_staged<T,I>::set_matrix_B(MatrixDenseView_new<T> * B_)
{
    if(B != nullptr) eslog::error("matrix B is already set\n");

    B = B_;
}



template<typename T, typename I>
void gemm_csx_dny_dny_staged<T,I>::set_matrix_C(MatrixDenseView_new<T> * C_)
{
    if(C != nullptr) eslog::error("matrix C is already set\n");

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
    stacktimer::push("gemm_csx_dny_dny_staged::preprocess");

    if(preprocess_called) eslog::error("preproces has already been called\n");
    if(A == nullptr) eslog::error("matrix A is not set\n");
    if(B == nullptr) eslog::error("matrix B is not set\n");
    if(C == nullptr) eslog::error("matrix C is not set\n");
    if(!A->ator->is_data_accessible_cpu()) eslog::error("matrix A must be cpu-accessible\n");
    if(!B->ator->is_data_accessible_cpu()) eslog::error("matrix B must be cpu-accessible\n");
    if(!C->ator->is_data_accessible_cpu()) eslog::error("matrix C must be cpu-accessible\n");
    if(A->nrows != C->nrows || B->ncols != C->ncols || A->ncols != B->nrows) eslog::error("incompatible matrix sizes\n");
    if(B->order != C->order) eslog::error("order of B and C must match\n");

    spblas::mm(*A, *B, *C, alpha, beta, handle_abc, 'P');

    stacktimer::pop();

    preprocess_called = true;
}



template<typename T, typename I>
void gemm_csx_dny_dny_staged<T,I>::perform()
{
    stacktimer::push("gemm_csx_dny_dny_staged::perform");

    if(!preprocess_called) eslog::error("preprocess has not been called\n");
    if(A == nullptr) eslog::error("matrix A is not set\n");
    if(B == nullptr) eslog::error("matrix B is not set\n");
    if(C == nullptr) eslog::error("matrix C is not set\n");
    if(A->nrows != C->nrows || B->ncols != C->ncols || A->ncols != B->nrows) eslog::error("incompatible matrix sizes\n");
    if(B->order != C->order) eslog::error("order of B and C must match\n");

    spblas::mm(*A, *B, *C, alpha, beta, handle_abc, 'C');

    stacktimer::pop();
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
        INSTANTIATE_T(std::complex<double>)

            INSTANTIATE

        #undef INSTANTIATE
    #undef INSTANTIATE_T
#undef INSTANTIATE_T_I



}
}
}
