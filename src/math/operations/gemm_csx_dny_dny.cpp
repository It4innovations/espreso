
#include "math/operations/gemm_csx_dny_dny.h"

#include "basis/utilities/stacktimer.h"



namespace espreso {
namespace math {
namespace operations {



template<typename T, typename I>
void gemm_csx_dny_dny<T,I>::set_matrix_A(MatrixCsxView_new<T,I> * A_)
{
    A = A_;
}



template<typename T, typename I>
void gemm_csx_dny_dny<T,I>::set_matrix_B(MatrixDenseView_new<T> * B_)
{
    B = B_;
}



template<typename T, typename I>
void gemm_csx_dny_dny<T,I>::set_matrix_C(MatrixDenseView_new<T> * C_)
{
    C = C_;
}



template<typename T, typename I>
void gemm_csx_dny_dny<T,I>::set_coefficients(T alpha_, T beta_)
{
    alpha = alpha_;
    beta = beta_;
}



template<typename T, typename I>
void gemm_csx_dny_dny<T,I>::perform()
{
    stacktimer::push("gemm_csx_dny_dny::perform");

    if(A == nullptr) eslog::error("matrix A is not set\n");
    if(B == nullptr) eslog::error("matrix B is not set\n");
    if(C == nullptr) eslog::error("matrix C is not set\n");
    if(A->nrows != C->nrows || B->ncols != C->ncols || A->ncols != B->nrows) eslog::error("incompatible matrix sizes\n");
    if(B->order != C->order) eslog::error("order of B and C must match\n");

    spblas::handle_mm handle_abc;
    spblas::mm(*A, *B, *C, alpha, beta, handle_abc, 'A');

    stacktimer::pop();
}



template<typename T, typename I>
void gemm_csx_dny_dny<T,I>::do_all(MatrixCsxView_new<T,I> * A, MatrixDenseView_new<T> * B, MatrixDenseView_new<T> * C, T alpha, T beta)
{
    gemm_csx_dny_dny<T,I> instance;
    instance.set_matrix_A(A);
    instance.set_matrix_B(B);
    instance.set_matrix_C(C);
    instance.set_coefficients(alpha, beta);
    instance.perform();
}



#define INSTANTIATE_T_I(T,I) \
template class gemm_csx_dny_dny<T,I>;

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
