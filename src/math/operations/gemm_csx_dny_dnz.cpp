
#include "math/operations/gemm_csx_dny_dnz.h"

#include "math/primitives_new/allocator_new.h"
#include "math/primitives_new/matrix_dense_data_new.h"
#include "math/operations/convert_dnx_dny.h"
#include "basis/utilities/stacktimer.h"



namespace espreso {
namespace math {
namespace operations {



template<typename T, typename I>
void gemm_csx_dny_dnz<T,I>::set_matrix_A(MatrixCsxView_new<T,I> * A_)
{
    if(A != nullptr) eslog::error("matrix A is already set\n");

    A = A_;
}



template<typename T, typename I>
void gemm_csx_dny_dnz<T,I>::set_matrix_B(MatrixDenseView_new<T> * B_)
{
    if(B != nullptr) eslog::error("matrix B is already set\n");

    B = B_;
}



template<typename T, typename I>
void gemm_csx_dny_dnz<T,I>::set_matrix_C(MatrixDenseView_new<T> * C_)
{
    if(C != nullptr) eslog::error("matrix C is already set\n");

    C = C_;
}



template<typename T, typename I>
void gemm_csx_dny_dnz<T,I>::set_coefficients(T alpha_, T beta_)
{
    alpha = alpha_;
    beta = beta_;
}



template<typename T, typename I>
void gemm_csx_dny_dnz<T,I>::perform()
{
    stacktimer::push("gemm_csx_dny_dnz::perform");

    if(A == nullptr) eslog::error("matrix A is not set\n");
    if(B == nullptr) eslog::error("matrix B is not set\n");
    if(C == nullptr) eslog::error("matrix C is not set\n");
    if(!A->ator->is_data_accessible_cpu()) eslog::error("matrix A must be cpu-accessible\n");
    if(!B->ator->is_data_accessible_cpu()) eslog::error("matrix B must be cpu-accessible\n");
    if(!C->ator->is_data_accessible_cpu()) eslog::error("matrix C must be cpu-accessible\n");
    if(A->nrows != C->nrows || B->ncols != C->ncols || A->ncols != B->nrows) eslog::error("incompatible matrix sizes\n");

    // spblas::mm requires B and C to have same order
    if(B->order == C->order) {
        spblas::handle_mm handle_abc;
        spblas::mm(*A, *B, *C, alpha, beta, handle_abc, 'A');
    }
    else {
        MatrixDenseData_new<T> C_tmp;
        C_tmp.set(C->nrows, C->ncols, change_order(C->order), AllocatorCPU_new::get_singleton());
        C_tmp.alloc();

        spblas::handle_mm handle_abc;
        spblas::mm(*A, *B, C_tmp, alpha, beta, handle_abc, 'A');

        convert_dnx_dny<T>::do_all(&C_tmp, C, false);
    }

    stacktimer::pop();
}



template<typename T, typename I>
void gemm_csx_dny_dnz<T,I>::do_all(MatrixCsxView_new<T,I> * A, MatrixDenseView_new<T> * B, MatrixDenseView_new<T> * C, T alpha, T beta)
{
    gemm_csx_dny_dnz<T,I> instance;
    instance.set_matrix_A(A);
    instance.set_matrix_B(B);
    instance.set_matrix_C(C);
    instance.set_coefficients(alpha, beta);
    instance.perform();
}



#define INSTANTIATE_T_I(T,I) \
template class gemm_csx_dny_dnz<T,I>;

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
