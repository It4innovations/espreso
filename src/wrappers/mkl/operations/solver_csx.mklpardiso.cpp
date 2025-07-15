
#ifdef HAVE_MKL

#include "wrappers/mkl/operations/solver_csx.mklpardiso.h"

#include "wrappers/mkl/pardiso_common.h"
#include "math/primitives_new/allocator_new.h"
#include "math/operations/convert_dnx_dny.h"



namespace espreso {
namespace math {
namespace operations {



template<typename T, typename I>
struct solver_csx_mklpardiso_data
{
    void * pt[64];
    MKL_INT iparm[64];
    MKL_INT mtype;
    MKL_INT msglvl;
    MKL_INT error;
    MKL_INT size_matrix_mklint;
};



template<typename T, typename I>
solver_csx_mklpardiso<T,I>::solver_csx_mklpardiso()
{
    if constexpr(!std::is_same_v<I,MKL_INT>) eslog::error("I and MKL_INT types dont match\n");

    data = std::make_unique<solver_csx_mklpardiso_data<T,I>>();
}



template<typename T, typename I>
solver_csx_mklpardiso<T,I>::~solver_csx_mklpardiso()
{
}



template<typename T, typename I>
void solver_csx_mklpardiso<T,I>::internal_factorize_symbolic()
{
    if(need_factors) eslog::error("no support for factor extraction\n");

    data->size_matrix_mklint = A->nrows;

    std::fill_n(data->pt, 64, nullptr);
    std::fill_n(data->iparm, 64, 0);
    data->iparm[0] = 1; // I did popullate iparm
    data->iparm[1] = 2; // metis fill reduction
    data->iparm[7] = 0; // Max number of refinement iterations (default changed between 2024.2 and 2025.0)
    data->iparm[34] = 1; // zero-based indexing

    data->mtype = get_pardiso_matrix_type<T>(A->prop);

    MKL_INT phase = 11;
    MKL_INT one = 1;
    MKL_INT nrhs = 0;
    pardiso(data->pt, &one, &one, &data->mtype, &phase, &data->size_matrix_mklint, A->vals, A->ptrs, A->idxs, nullptr, &nrhs, data->iparm, &data->msglvl, nullptr, nullptr, &data->error);
    if(data->error != 0) {
        eslog::error("pardiso error %d: %s\n", (int)data->error, get_pardiso_error_string(data->error));
    }
}



template<typename T, typename I>
void solver_csx_mklpardiso<T,I>::internal_factorize_numeric()
{
    MKL_INT phase = 22;
    MKL_INT one = 1;
    MKL_INT nrhs = 0;
    pardiso(data->pt, &one, &one, &data->mtype, &phase, &data->size_matrix_mklint, A->vals, A->ptrs, A->idxs, nullptr, &nrhs, data->iparm, &data->msglvl, nullptr, nullptr, &data->error);
    if(data->error != 0) {
        eslog::error("pardiso error %d: %s\n", (int)data->error, get_pardiso_error_string(data->error));
    }
}



template<typename T, typename I>
void solver_csx_mklpardiso<T,I>::internal_get_permutation(PermutationView_new<I> & perm)
{
    eslog::error("no support for factor extraction\n");
}



template<typename T, typename I>
void solver_csx_mklpardiso<T,I>::internal_get_factor_L(MatrixCsxView_new<T,I> & L, bool pattern, bool values)
{
    eslog::error("no support for factor extraction\n");
}



template<typename T, typename I>
void solver_csx_mklpardiso<T,I>::internal_get_factor_U(MatrixCsxView_new<T,I> & U, bool pattern, bool values)
{
    eslog::error("no support for factor extraction\n");
}



template<typename T, typename I>
void solver_csx_mklpardiso<T,I>::internal_solve(VectorDenseView_new<T> & rhs, VectorDenseView_new<T> & sol)
{
    MKL_INT phase = 33;
    MKL_INT one = 1;
    MKL_INT nrhs = 1;
    pardiso(data->pt, &one, &one, &data->mtype, &phase, &data->size_matrix_mklint, A->vals, A->ptrs, A->idxs, nullptr, &nrhs, data->iparm, &data->msglvl, rhs.vals, sol.vals, &data->error);
    if(data->error != 0) {
        eslog::error("pardiso error %d: %s\n", (int)data->error, get_pardiso_error_string(data->error));
    }
}



template<typename T, typename I>
void solver_csx_mklpardiso<T,I>::internal_solve(MatrixDenseView_new<T> & rhs, MatrixDenseView_new<T> & sol)
{
    // no support for leading dimension in pardiso. matrix has to be colmajor. so I have to reallocate and convert
    Allocator_new * ator = AllocatorCPU_new::get_singleton();
    T * rhs_nold_vals = ator->template alloc<T>(rhs.nrows * rhs.ncols);
    T * sol_nold_vals = ator->template alloc<T>(sol.nrows * sol.ncols);

    {
        MatrixDenseView_new<T> rhs_nold;
        rhs_nold.set_view(rhs.nrows, rhs.ncols, rhs.nrows, 'C', rhs_nold_vals, ator);

        MatrixDenseView_new<T> sol_nold;
        sol_nold.set_view(sol.nrows, sol.ncols, sol.nrows, 'C', sol_nold_vals, ator);

        convert_dnx_dny<T>::do_all(&rhs, &rhs_nold, false);

        MKL_INT phase = 33;
        MKL_INT one = 1;
        MKL_INT nrhs = rhs.ncols;
        pardiso(data->pt, &one, &one, &data->mtype, &phase, &data->size_matrix_mklint, A->vals, A->ptrs, A->idxs, nullptr, &nrhs, data->iparm, &data->msglvl, rhs_nold.vals, sol_nold.vals, &data->error);
        if(data->error != 0) {
            eslog::error("pardiso error %d: %s\n", (int)data->error, get_pardiso_error_string(data->error));
        }
        
        convert_dnx_dny<T>::do_all(&sol_nold, &sol, false);
    }

    ator->free(rhs_nold_vals);
    ator->free(sol_nold_vals);
}



#define INSTANTIATE_T_I(T,I) \
template class solver_csx_mklpardiso<T,I>;

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

#endif
