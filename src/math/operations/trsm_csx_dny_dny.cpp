
#include "math/operations/trsm_csx_dny_dny.h"

#include "math/wrappers/math.spblas.h"
#include "math/primitives_new/allocator_new.h"
#include "math/operations/copy_dnx.h"
#include "basis/utilities/stacktimer.h"



namespace espreso {
namespace math {
namespace operations {



template<typename T, typename I>
void trsm_csx_dny_dny<T,I>::set_system_matrix(MatrixCsxView_new<T,I> * A_)
{
    if(A != nullptr) eslog::error("matrix A is already set\n");

    A = A_;
}



template<typename T, typename I>
void trsm_csx_dny_dny<T,I>::set_rhs_matrix(MatrixDenseView_new<T> * B_)
{
    if(B != nullptr) eslog::error("matrix B is already set\n");

    B = B_;
}



template<typename T, typename I>
void trsm_csx_dny_dny<T,I>::set_solution_matrix(MatrixDenseView_new<T> * X_)
{
    if(X != nullptr) eslog::error("matrix X is already set\n");

    X = X_;
}



template<typename T, typename I>
void trsm_csx_dny_dny<T,I>::perform()
{
    stacktimer::push("trsm_csx_dny_dny::perform");

    if(A == nullptr) eslog::error("system matrix is not set\n");
    if(X == nullptr) eslog::error("solution matrix is not set\n");
    if(B == nullptr) eslog::error("rhs matrix is not set\n");
    if(!A->ator->is_data_accessible_cpu()) eslog::error("matrix A must be cpu-accessible\n");
    if(!B->ator->is_data_accessible_cpu()) eslog::error("matrix B must be cpu-accessible\n");
    if(!X->ator->is_data_accessible_cpu()) eslog::error("matrix X must be cpu-accessible\n");
    if(A->nrows != A->ncols) eslog::error("system matrix must be square\n");
    if(X->nrows != B->nrows || X->ncols != B->ncols) eslog::error("rhs and sol matrix sizes dont match\n");
    if(X->order != B->order) eslog::error("rhs and sol orders dont match\n");
    if(A->nrows != X->nrows) eslog::error("incompatible matrices\n");
    if(A->prop.uplo != 'L' && A->prop.uplo != 'U') eslog::error("matrix A has wrong uplo\n");

    if(X == B) { // in-place
        MatrixDenseData_new<T> Y;
        Y.set(X->nrows, X->ncols, X->order, AllocatorCPU_new::get_singleton());
        Y.alloc();
        trsm_csx_dny_dny<T,I>::do_all(A, &Y, B);
        copy_dnx<T>::do_all(&Y, X);
        Y.clear();
    }
    else {
        spblas::handle_trsm handle;
        math::spblas::trsm(*A, *X, *B, handle, 'A');
    }

    stacktimer::pop();
}



template<typename T, typename I>
void trsm_csx_dny_dny<T,I>::do_all(MatrixCsxView_new<T,I> * A, MatrixDenseView_new<T> * B, MatrixDenseView_new<T> * X)
{
    trsm_csx_dny_dny<T,I> instance;
    instance.set_system_matrix(A);
    instance.set_rhs_matrix(B);
    instance.set_solution_matrix(X);
    instance.perform();
}



#define INSTANTIATE_T_I(T,I) \
template class trsm_csx_dny_dny<T,I>;

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
