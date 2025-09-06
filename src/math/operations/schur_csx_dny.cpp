
#include "math/operations/schur_csx_dny.h"

#include "math/operations/schur_csx_dny.tria.h"
#include "math/operations/schur_csx_dny.spsolver.h"
#include "wrappers/mkl/operations/schur_csx_dny.mklpardiso.h"
#include "wrappers/mumps/operations/schur_csx_dny.mumps.h"
#include "wrappers/pastix/operations/schur_csx_dny.pastix.h"
#include "math/primitives_new/allocator_new.h"
#include "math/operations/concat_csx.h"
#include "basis/utilities/stacktimer.h"



namespace espreso {
namespace math {
namespace operations {



template<typename T, typename I>
std::unique_ptr<schur_csx_dny<T,I>> schur_csx_dny<T,I>::make(implementation_selector is)
{
    auto autoselect_implementation = [](){
        #ifdef HAVE_PASTIX
            return implementation_selector::pastix;
        #endif
        #ifdef HAVE_MKL
            return implementation_selector::mklpardiso;
        #endif
        return implementation_selector::triangular;
    };

    switch(is) {
        case implementation_selector::autoselect:
            return schur_csx_dny<T,I>::make(autoselect_implementation());
        case implementation_selector::triangular:
            return std::make_unique<schur_csx_dny_tria<T,I>>();
        case implementation_selector::mklpardiso:
            return std::make_unique<schur_csx_dny_mklpardiso<T,I>>();
        case implementation_selector::sparse_solver:
            return std::make_unique<schur_csx_dny_spsolver<T,I>>();
        case implementation_selector::mumps:
            return std::make_unique<schur_csx_dny_mumps<T,I>>();
        case implementation_selector::pastix:
            return std::make_unique<schur_csx_dny_pastix<T,I>>();
        default:
            eslog::error("invalid implementation selector\n");
    }
}



template<typename T, typename I>
void schur_csx_dny<T,I>::set_coefficients(Treal alpha_)
{
    alpha = alpha_;
}



template<typename T, typename I>
void schur_csx_dny<T,I>::set_matrix(MatrixCsxView_new<T,I> * A11_, MatrixCsxView_new<T,I> * A12_, MatrixCsxView_new<T,I> * A21_, MatrixCsxView_new<T,I> * A22_)
{
    if(called_set_matrix != '_') eslog::error("matrix is already set\n");

    A11 = A11_;
    A12 = A12_;
    A21 = A21_;
    A22 = A22_;

    called_set_matrix = '4';
}



template<typename T, typename I>
void schur_csx_dny<T,I>::set_matrix(MatrixCsxView_new<T,I> * A_, size_t size_sc_)
{
    if(called_set_matrix != '_') eslog::error("matrix is already set\n");

    A = A_;
    size_sc = size_sc_;

    called_set_matrix = '1';
}



template<typename T, typename I>
void schur_csx_dny<T,I>::set_sc(MatrixDenseView_new<T> * sc_)
{
    if(sc != nullptr) eslog::error("matrix sc is already set\n");

    sc = sc_;
}



template<typename T, typename I>
void schur_csx_dny<T,I>::set_need_solve_A11(bool need_solve_A11_)
{
    if(called_preprocess) eslog::error("cannot re-set need_solve_A11 after preprocess was called\n");

    need_solve_A11 = need_solve_A11_;
}



template<typename T, typename I>
void schur_csx_dny<T,I>::preprocess()
{
    stacktimer::push("schur_csx_dny::preprocess");

    if(called_preprocess) eslog::error("preprocess was already called\n");
    if(called_set_matrix == '_') eslog::error("matrix is not set\n");
    if(sc == nullptr) eslog::error("sc is not set\n");
    if(!sc->ator->is_data_accessible_cpu()) eslog::error("matrix sc must be cpu-accessible\n");
    if(sc->nrows != sc->ncols) eslog::error("sc has to be square\n");

    if(called_set_matrix == '4') {
        if(A11 != nullptr && !A11->ator->is_data_accessible_cpu()) eslog::error("matrix A11 must be cpu-accessible\n");
        if(A12 != nullptr && !A12->ator->is_data_accessible_cpu()) eslog::error("matrix A12 must be cpu-accessible\n");
        if(A21 != nullptr && !A21->ator->is_data_accessible_cpu()) eslog::error("matrix A21 must be cpu-accessible\n");
        if(A22 != nullptr && !A22->ator->is_data_accessible_cpu()) eslog::error("matrix A22 must be cpu-accessible\n");
        if(A11 == nullptr) eslog::error("A11 cannot be nullptr\n");
        int num_sides_set = (int)(A12 != nullptr) + (int)(A21 != nullptr);
        if(is_hermitian<T>(A11->prop.symm) && num_sides_set == 0) eslog::error("at least one of A12 and A21 has to be set\n");
        if(!is_hermitian<T>(A11->prop.symm) && num_sides_set <= 1) eslog::error("both A12 and A21 have to be set\n");
        if(A11->nrows != A11->ncols) eslog::error("A11 has to be square\n");
        if(A22 != nullptr && A22->nrows != A22->ncols) eslog::error("A22 has to be square\n");
        if(A22 != nullptr && A22->prop.symm != A11->prop.symm) eslog::error("A11 and A22 must have equal symmetry\n");
        if(A11->prop.symm != sc->prop.symm) eslog::error("A11 and sc must have equal symmetry\n");
        if(A12 != nullptr && A12->nrows != A11->nrows) eslog::error("wrong matrix A12 size (does not match A11)\n");
        if(A21 != nullptr && A21->ncols != A11->ncols) eslog::error("wrong matrix A21 size (does not match A11)\n");
        if(A12 != nullptr && A22 != nullptr && A12->ncols != A22->ncols) eslog::error("wrong matrix A12 size (does not natch A22)\n");
        if(A21 != nullptr && A22 != nullptr && A21->nrows != A22->nrows) eslog::error("wrong matrix A21 size (does not natch A22)\n");
        // if only one of A12 or A21 is set and the system is hermitian, it is assumed the other is its conjugate transpose
        // if A22 is nullptr, it is considered to be a null matrix

        if(A12 != nullptr) size_sc = A12->ncols;
        if(A21 != nullptr) size_sc = A21->nrows;
        size_A11 = A11->nrows;
        size_matrix = size_A11 + size_sc;

        is_matrix_hermitian = is_hermitian<T>(A11->prop.symm);

        if(is_matrix_hermitian && A11->prop.uplo != 'U' && A11->prop.uplo != 'L') eslog::error("wrong A11 uplo\n");
        if(is_matrix_hermitian && sc->prop.uplo != 'U' && sc->prop.uplo != 'L') eslog::error("wrong sc uplo\n");
        if(is_matrix_hermitian && A22 != nullptr && A22->prop.uplo != 'U' && A22->prop.uplo != 'L') eslog::error("wrong A22 uplo\n");
    }

    if(called_set_matrix == '1') {
        if(!A->ator->is_data_accessible_cpu()) eslog::error("matrix A must be cpu-accessible\n");

        size_matrix = A->nrows;
        size_A11 = size_matrix - size_sc;

        if(is_symmetric<T>(A->prop.symm) != is_symmetric<T>(sc->prop.symm) || is_hermitian<T>(A->prop.symm) != is_hermitian<T>(sc->prop.symm)) eslog::error("A and sc must have equivalent symmetry\n");

        is_matrix_hermitian = is_hermitian<T>(A->prop.symm);
        
        if(is_matrix_hermitian && A->prop.uplo != 'U' && A->prop.uplo != 'L') eslog::error("wrong A uplo\n");
        if(is_matrix_hermitian && sc->prop.uplo != 'U' && sc->prop.uplo != 'L') eslog::error("wrong sc uplo\n");
    }

    if(sc->nrows != size_sc) eslog::error("mismatch between provided size_sc and SC matrix size\n");

    this->internal_preprocess();

    called_preprocess = true;

    stacktimer::pop();
}



template<typename T, typename I>
void schur_csx_dny<T,I>::perform_1()
{
    stacktimer::push("schur_csx_dny::perform_1");

    if(!called_preprocess) eslog::error("preprocess has not been called\n");

    this->internal_perform_1();

    stacktimer::pop();
}



template<typename T, typename I>
void schur_csx_dny<T,I>::perform_2()
{
    stacktimer::push("schur_csx_dny::perform_2");

    if(!called_preprocess) eslog::error("preprocess has not been called\n");

    this->internal_perform_2();

    called_perform = true;

    stacktimer::pop();
}



template<typename T, typename I>
void schur_csx_dny<T,I>::perform()
{
    if(!called_preprocess) eslog::error("preprocess has not been called\n");

    perform_1();
    perform_2();
}



template<typename T, typename I>
void schur_csx_dny<T,I>::solve_A11(VectorDenseView_new<T> & rhs, VectorDenseView_new<T> & sol)
{
    stacktimer::push("schur_csx_dny::solve_A11");

    if(!called_preprocess) eslog::error("preprocess has not been called\n");
    if(!need_solve_A11) eslog::error("need_solve_A11 is not set, so cannot solve A11\n");
    if(rhs.size != size_A11) eslog::error("wrong rhs size\n");
    if(sol.size != size_A11) eslog::error("wrong sol size\n");

    this->internal_solve_A11(rhs, sol);

    stacktimer::pop();
}



template<typename T, typename I>
void schur_csx_dny<T,I>::helper_concat(MatrixCsxData_new<T,I> & A_whole, char stage)
{
    if(called_set_matrix != '4') eslog::error("it does not make sense to concat\n");

    if(stage == 'P') {
        if(A != nullptr) eslog::error("helper_concat preprocess already caller, or other error\n");

        if(is_matrix_hermitian) {
            if(A11->prop.uplo != 'U') eslog::error("only uplo=U is supported\n");
            if(A22 != nullptr && A22->prop.uplo != 'U') eslog::error("only uplo=U is supported\n");
            if(A21 != nullptr) eslog::error("A21 should be empty\n");
            if(A11->prop.dfnt != MatrixDefinitness_new::positive_definite) eslog::error("A11 should be positive definite\n");
        }
        size_t total_nnz = 0;
        if(A11 != nullptr) total_nnz += A11->nnz;
        if(A12 != nullptr) total_nnz += A12->nnz;
        if(A21 != nullptr) total_nnz += A21->nnz;
        if(A22 != nullptr) total_nnz += A22->nnz;
        A_whole.set(size_matrix, size_matrix, total_nnz, 'R', AllocatorCPU_new::get_singleton());
        A_whole.alloc();
        if(is_matrix_hermitian) {
            A_whole.prop.uplo = 'U';
            A_whole.prop.symm = MatrixSymmetry_new::hermitian;
            A_whole.prop.dfnt = MatrixDefinitness_new::positive_semidefinite;
        }
        else {
            A_whole.prop.uplo = 'F';
            A_whole.prop.symm = MatrixSymmetry_new::general;
            A_whole.prop.dfnt = MatrixDefinitness_new::indefinite;
        }
        A = &A_whole;
        concat_csx<T,I>::do_all(A11, A12, A21, A22, A);
    }

    if(stage == 'F') {
        if(A == nullptr) eslog::error("helper_concat preprocess not called\n");
        concat_csx<T,I>::do_all(A11, A12, A21, A22, A);
    }
}



#define INSTANTIATE_T_I(T,I) \
template class schur_csx_dny<T,I>;

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
