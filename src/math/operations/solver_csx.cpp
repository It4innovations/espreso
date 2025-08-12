
#include "math/operations/solver_csx.h"

#include "wrappers/suitesparse/operations/solver_csx.cholmod.h"
#include "wrappers/suitesparse/operations/solver_csx.umfpack.h"
#include "wrappers/mkl/operations/solver_csx.mklpardiso.h"
#include "wrappers/mumps/operations/solver_csx.mumps.h"
#include "wrappers/strumpack/operations/solver_csx.strumpack.h"
#include "wrappers/pastix/operations/solver_csx.pastix.h"
#include "basis/utilities/stacktimer.h"



namespace espreso {
namespace math {
namespace operations {



template<typename T, typename I>
std::unique_ptr<solver_csx<T,I>> solver_csx<T,I>::make(implementation_selector is, MatrixCsxView_new<T,I> * matrix, bool need_factors, bool need_solve)
{
    std::vector<implementation_selector> implementations;
    if(is == implementation_selector::autoselect) {
        #ifdef HAVE_MKL
            implementations.push_back(implementation_selector::mklpardiso);
        #endif
        #ifdef HAVE_SUITESPARSE
            implementations.push_back(implementation_selector::suitesparse);
        #endif
        #ifdef HAVE_STRUMPACK
            implementations.push_back(implementation_selector::strumpack);
        #endif
        #ifdef HAVE_MUMPS
            implementations.push_back(implementation_selector::mumps);
        #endif
        #ifdef HAVE_PASTIX
            implementations.push_back(implementation_selector::pastix);
        #endif
    }
    else {
        implementations.push_back(is);
    }

    for(implementation_selector impl : implementations) {
        switch(impl) {
            case implementation_selector::mklpardiso:
                if(!need_factors) return std::make_unique<solver_csx_mklpardiso<T,I>>();
                break;
            case implementation_selector::suitesparse:
                if(matrix != nullptr && is_hermitian<T>(matrix->prop.symm)) return std::make_unique<solver_csx_cholmod<T,I>>();
                else return std::make_unique<solver_csx_umfpack<T,I>>();
                break;
            case implementation_selector::mumps:
                if(!need_factors) return std::make_unique<solver_csx_mumps<T,I>>();
                break;
            case implementation_selector::strumpack:
                if(!need_factors) return std::make_unique<solver_csx_strumpack<T,I>>();
                break;
            case implementation_selector::pastix:
                if(!need_factors) return std::make_unique<solver_csx_pastix<T,I>>();
                break;
            default:
                eslog::error("invalid solver_csx implementation selector\n");
        }
    }

    if(is == implementation_selector::autoselect) {
        eslog::error("no valid solver_csx implementation available\n");
    }
    eslog::error("requested implementation does not support requested matrix and needs\n");
}



template<typename T, typename I>
void solver_csx<T,I>::set_matrix_A(MatrixCsxView_new<T,I> * A_)
{
    if(A != nullptr) eslog::error("matrix A is already set\n");

    A = A_;
}



template<typename T, typename I>
void solver_csx<T,I>::set_needs(bool need_factors_, bool need_solve_)
{
    if(called_factorize_symbolic) eslog::error("cannot change needs after factorization\n");

    need_factors = need_factors_;
    need_solve = need_solve_;
}



template<typename T, typename I>
void solver_csx<T,I>::factorize_symbolic()
{
    stacktimer::push("solver_csx::factorize_symbolic");

    if(called_factorize_symbolic) eslog::error("symbolic factorization has already been called\n");
    if(A == nullptr) eslog::error("matrix A is not set\n");
    if(A->nrows != A->ncols) eslog::error("A has to be square\n");
    if(is_hermitian<T>(A->prop.symm) && A->prop.uplo != 'L' && A->prop.uplo != 'U') eslog::error("unset A uplo\n");
    if(!need_factors && !need_solve) eslog::error("nonsensical needs\n");

    this->internal_factorize_symbolic();

    called_factorize_symbolic = true;

    stacktimer::pop();
}



template<typename T, typename I>
void solver_csx<T,I>::factorize_numeric()
{
    stacktimer::push("solver_csx::factorize_numeric");

    if(!called_factorize_symbolic) eslog::error("symbolic factorization has not been called\n");

    this->internal_factorize_numeric();

    called_factorize_numeric = true;

    stacktimer::pop();
}



template<typename T, typename I>
size_t solver_csx<T,I>::get_factor_nnz_L()
{
    if(!need_factors) eslog::error("need factors was not set\n");
    if(!called_factorize_symbolic) eslog::error("symbolic factorization has not been called\n");

    return nnz_L;
}



template<typename T, typename I>
size_t solver_csx<T,I>::get_factor_nnz_U()
{
    if(!need_factors) eslog::error("need factors was not set\n");
    if(!called_factorize_symbolic) eslog::error("symbolic factorization has not been called\n");

    return nnz_U;
}



template<typename T, typename I>
size_t solver_csx<T,I>::get_factor_nnz(char uplo)
{
    if(!called_factorize_symbolic) eslog::error("symbolic factorization has not been called\n");

    if(uplo == 'L') return nnz_L;
    if(uplo == 'U') return nnz_U;
    if(A->prop.uplo == 'L') return nnz_L;
    if(A->prop.uplo == 'U') return nnz_U;
    eslog::error("wrong call to get_factor_nnz\n");
}



template<typename T, typename I>
void solver_csx<T,I>::get_permutation(PermutationView_new<I> & perm)
{
    if(!called_factorize_symbolic) eslog::error("numeric factorization has not been called\n");
    if(perm.size != A->nrows) eslog::error("wrong permutation size\n");

    this->internal_get_permutation(perm);
}



template<typename T, typename I>
void solver_csx<T,I>::get_factor_L(MatrixCsxView_new<T,I> & L, bool pattern, bool values)
{
    stacktimer::push("solver_csx::get_factor_L");

    if(!need_factors) eslog::error("need factors was not set\n");
    if(!called_factorize_symbolic) eslog::error("symbolic factorization has not been called\n");
    if(values && !called_factorize_numeric) eslog::error("numeric factorization has not been called\n");
    if(L.order != A->order) eslog::error("order of L and A must match\n");
    if(L.nrows != A->nrows || L.ncols != A->ncols || L.nnz != nnz_L) eslog::error("wrong L matrix dimensions\n");
    if(L.prop.uplo != 'L') eslog::error("wrong factor uplo\n");
    if(is_hermitian<T>(A->prop.symm) && A->prop.uplo != 'L') eslog::error("querrying incorrect factor\n");

    this->internal_get_factor_L(L, pattern, values);

    stacktimer::pop();
}



template<typename T, typename I>
void solver_csx<T,I>::get_factor_U(MatrixCsxView_new<T,I> & U, bool pattern, bool values)
{
    stacktimer::push("solver_csx::get_factor_U");

    if(!need_factors) eslog::error("need factors was not set\n");
    if(!called_factorize_symbolic) eslog::error("symbolic factorization has not been called\n");
    if(values && !called_factorize_numeric) eslog::error("numeric factorization has not been called\n");
    if(U.order != A->order) eslog::error("order of U and A must match\n");
    if(U.nrows != A->nrows || U.ncols != A->ncols || U.nnz != nnz_U) eslog::error("wrong U matrix dimensions\n");
    if(U.prop.uplo != 'U') eslog::error("wrong factor uplo\n");
    if(is_hermitian<T>(A->prop.symm) && A->prop.uplo != 'U') eslog::error("querrying incorrect factor\n");

    this->internal_get_factor_U(U, pattern, values);

    stacktimer::pop();
}



template<typename T, typename I>
void solver_csx<T,I>::get_factor(MatrixCsxView_new<T,I> & factor, bool pattern, bool values)
{
    if(factor.prop.uplo == 'L') {
        get_factor_L(factor, pattern, values);
    }
    else if(factor.prop.uplo == 'U') {
        get_factor_U(factor, pattern, values);
    }
    else {
        eslog::error("unset factor uplo\n");
    }
}



template<typename T, typename I>
void solver_csx<T,I>::solve(VectorDenseView_new<T> & rhs, VectorDenseView_new<T> & sol)
{
    stacktimer::push("solver_csx::solve (vector)");

    if(!need_solve) eslog::error("need solve was not set\n");
    if(!called_factorize_numeric) eslog::error("numeric factorization has not been called\n");
    if(rhs.size != A->nrows) eslog::error("wrong rhs size\n");
    if(sol.size != A->nrows) eslog::error("wrong sol size\n");

    this->internal_solve(rhs, sol);

    stacktimer::pop();
}



template<typename T, typename I>
void solver_csx<T,I>::solve(MatrixDenseView_new<T> & rhs, MatrixDenseView_new<T> & sol)
{
    stacktimer::push("solver_csx::solve (matrix dense)");

    if(!need_solve) eslog::error("need solve was not set\n");
    if(!called_factorize_numeric) eslog::error("numeric factorization has not been called\n");
    if(rhs.nrows != A->nrows) eslog::error("wrong rhs nrows\n");
    if(sol.nrows != A->nrows) eslog::error("wrong sol nrows\n");
    if(rhs.ncols != sol.ncols) eslog::error("rhs and sol ncols does not match\n");
    if(rhs.order != sol.order) eslog::error("rhs and sol orders dont match\n");

    this->internal_solve(rhs, sol);

    stacktimer::pop();
}



template<typename T, typename I>
void solver_csx<T,I>::solve(MatrixCsxView_new<T,I> & rhs, MatrixDenseView_new<T> & sol)
{
    stacktimer::push("solver_csx::solve (matrix sparse)");

    if(!need_solve) eslog::error("need solve was not set\n");
    if(!called_factorize_numeric) eslog::error("numeric factorization has not been called\n");
    if(rhs.nrows != A->nrows) eslog::error("wrong rhs nrows\n");
    if(sol.nrows != A->nrows) eslog::error("wrong sol nrows\n");
    if(rhs.ncols != sol.ncols) eslog::error("rhs and sol ncols does not match\n");

    this->internal_solve(rhs, sol);

    stacktimer::pop();
}



#define INSTANTIATE_T_I(T,I) \
template class solver_csx<T,I>;

    #define INSTANTIATE_T(T) \
    INSTANTIATE_T_I(T, int32_t) \
    /* INSTANTIATE_T_I(T, int64_t) */

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
