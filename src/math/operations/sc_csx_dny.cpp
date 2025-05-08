
#include "math/operations/sc_csx_dny.h"

#include "math/operations/sc_csx_dny.tria.h"
#include "math/operations/sc_csx_dny.spsolver.h"
#include "wrappers/mkl/operations/sc_csx_dny.mklpardiso.h"
#include "basis/utilities/stacktimer.h"



namespace espreso {
namespace math {
namespace operations {



template<typename T, typename I>
std::unique_ptr<sc_csx_dny<T,I>> sc_csx_dny<T,I>::make(implementation_selector is)
{
    auto autoselect_implementation = [](){
        #ifdef HAVE_MKL
            return implementation_selector::mklpardiso;
        #endif
        return implementation_selector::triangular;
    };

    switch(is) {
        case implementation_selector::autoselect:
            return sc_csx_dny<T,I>::make(autoselect_implementation());
        case implementation_selector::triangular:
            return std::make_unique<sc_csx_dny_tria<T,I>>();
        case implementation_selector::mklpardiso:
            return std::make_unique<sc_csx_dny_mklpardiso<T,I>>();
        case implementation_selector::spsolver:
            return std::make_unique<sc_csx_dny_spsolver<T,I>>();
        default:
            eslog::error("invalid implementation selector\n");
    }
}



template<typename T, typename I>
void sc_csx_dny<T,I>::set_coefficients(Treal alpha_)
{
    alpha = alpha_;
}



template<typename T, typename I>
void sc_csx_dny<T,I>::set_matrix(MatrixCsxView_new<T,I> * A11_, MatrixCsxView_new<T,I> * A12_, MatrixCsxView_new<T,I> * A21_, MatrixCsxView_new<T,I> * A22_)
{
    if(called_set_matrix != '_') eslog::error("matrix is already set\n");

    A11 = A11_;
    A12 = A12_;
    A21 = A21_;
    A22 = A22_;

    called_set_matrix = '4';
}



template<typename T, typename I>
void sc_csx_dny<T,I>::set_matrix(MatrixCsxView_new<T,I> * A_, size_t size_sc_)
{
    if(called_set_matrix != '_') eslog::error("matrix is already set\n");

    A = A_;
    size_sc = size_sc_;

    called_set_matrix = '1';
}



template<typename T, typename I>
void sc_csx_dny<T,I>::set_sc(MatrixDenseView_new<T> * sc_)
{
    if(sc != nullptr) eslog::error("matrix sc is already set\n");
    if(sc_ == nullptr) eslog::error("sc cannot be nullptr\n");

    sc = sc_;

    if(sc->nrows != sc->ncols) eslog::error("sc has to be square\n");
}



template<typename T, typename I>
void sc_csx_dny<T,I>::set_need_solve_A11(bool need_solve_A11_)
{
    need_solve_A11 = need_solve_A11_;
}



template<typename T, typename I>
void sc_csx_dny<T,I>::preprocess()
{
    stacktimer::push("sc_csx_dny::preprocess");

    if(called_preprocess) eslog::error("preprocess was already called\n");
    if(called_set_matrix == '_') eslog::error("matrix is not set\n");
    if(sc == nullptr) eslog::error("sc is not set\n");

    if(called_set_matrix == '4') {
        if(A11 == nullptr) eslog::error("A11 cannot be nullptr\n");
        int num_sides_set = (int)(A12 != nullptr) + (int)(A21 != nullptr);
        if(is_hermitian<T>(A11->prop.symm) && num_sides_set == 0) eslog::error("at least one of A12 and A21 has to be set\n");
        if(!is_hermitian<T>(A11->prop.symm) && num_sides_set <= 1) eslog::error("both A12 and A21 have to be set\n");
        if(A11->nrows != A11->ncols) eslog::error("A11 has to be square\n");
        if(A22 != nullptr && A22->nrows != A22->ncols) eslog::error("A22 has to be square\n");
        if(A22 != nullptr && A22->prop.symm != A11->prop.symm) eslog::error("A11 and A22 must have equal symmetry\n");
        if(A11->prop.symm != sc->prop.symm) eslog::error("A11 and sc must have equal symmetry\n");
        if(A12 != nullptr && A12->nrows != A11->nrows) eslog::error("wrong matrix A12 size (does not natch A11)\n");
        if(A21 != nullptr && A21->ncols != A11->ncols) eslog::error("wrong matrix A21 size (does not natch A11)\n");
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
        size_matrix = A->nrows;
        size_A11 = size_matrix - size_sc;

        if(A->prop.symm != sc->prop.symm) eslog::error("A and sc must have equal symmetry\n");

        is_matrix_hermitian = is_hermitian<T>(A11->prop.symm);
        
        if(is_matrix_hermitian && A->prop.uplo != 'U' && A->prop.uplo != 'L') eslog::error("wrong A uplo\n");
        if(is_matrix_hermitian && sc->prop.uplo != 'U' && sc->prop.uplo != 'L') eslog::error("wrong sc uplo\n");
    }

    this->internal_preprocess();

    called_preprocess = true;

    stacktimer::pop();
}



template<typename T, typename I>
void sc_csx_dny<T,I>::perform_1()
{
    stacktimer::push("sc_csx_dny::perform_1");

    if(!called_preprocess) eslog::error("preprocess has not been called\n");

    this->internal_perform_1();

    stacktimer::pop();
}



template<typename T, typename I>
void sc_csx_dny<T,I>::perform_2()
{
    stacktimer::push("sc_csx_dny::perform_2");

    if(!called_preprocess) eslog::error("preprocess has not been called\n");

    this->internal_perform_2();

    called_perform = true;

    stacktimer::pop();
}



template<typename T, typename I>
void sc_csx_dny<T,I>::perform()
{
    if(!called_preprocess) eslog::error("preprocess has not been called\n");

    perform_1();
    perform_2();
}



template<typename T, typename I>
void sc_csx_dny<T,I>::solve_A11(VectorDenseView_new<T> & rhs, VectorDenseView_new<T> & sol)
{
    stacktimer::push("sc_csx_dny::solve_A11");

    if(!called_preprocess) eslog::error("preprocess has not been called\n");
    if(!need_solve_A11) eslog::error("need_solve_A11 is not set, so cannot solve A11\n");
    if(rhs.size != size_A11) eslog::error("wrong rhs size\n");
    if(sol.size != size_A11) eslog::error("wrong sol size\n");

    this->internal_solve_A11(rhs, sol);

    stacktimer::pop();
}



#define INSTANTIATE_T_I(T,I) \
template class sc_csx_dny<T,I>;

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
