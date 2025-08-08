
#ifdef HAVE_STRUMPACK

#include "wrappers/strumpack/operations/solver_csx.strumpack.h"

#include "wrappers/strumpack/strumpack_common.h"

#include "math/primitives_new/allocator_new.h"
#include "math/operations/convert_csx_dny.h"
#include "math/operations/complete_csx_csx_map.h"



namespace espreso {
namespace math {
namespace operations {



template<typename T, typename I>
struct solver_csx_strumpack_data
{
    strumpack::SparseSolver<T,I> solver = strumpack::SparseSolver<T,I>(false);
    complete_csx_csx_map<T,I> op_A_complete;
    MatrixCsxData_new<T,I> A_full;
    MatrixCsxView_new<T,I> * A_to_use;
    bool need_complete = false;
    bool called_factorize_numeric = false;
};



template<typename T, typename I>
solver_csx_strumpack<T,I>::solver_csx_strumpack()
{
    data = std::make_unique<solver_csx_strumpack_data<T,I>>();
}



template<typename T, typename I>
solver_csx_strumpack<T,I>::~solver_csx_strumpack()
{
}



template<typename T, typename I>
void solver_csx_strumpack<T,I>::internal_factorize_symbolic()
{
    if(A->order != 'R') eslog::error("support only CSR matrices for now\n");

    data->need_complete = is_uplo(A->prop.uplo);
    if(data->need_complete) {
        data->A_full.set(A->ncols, A->nrows, 2 * A->nnz - A->nrows, A->order, AllocatorCPU_new::get_singleton());
        data->A_full.alloc();
        data->A_full.prop = A->prop;
        data->A_full.prop.uplo = 'F';
        
        data->op_A_complete.set_matrix_src(A);
        data->op_A_complete.set_matrix_dst(&data->A_full);
        data->op_A_complete.set_conj(is_hermitian<T>(A->prop.symm));
        data->op_A_complete.perform_pattern();

        data->A_to_use = &data->A_full;
    }
    else {
        data->A_to_use = A;
    }

    if(is_symmetric<T>(data->A_to_use->prop.symm)) {
        data->solver.options().enable_symmetric();
    }
    if(data->A_to_use->prop.dfnt == MatrixDefinitness_new::positive_definite) {
        data->solver.options().enable_positive_definite();
    }
    data->solver.options().set_Krylov_solver(strumpack::KrylovSolver::DIRECT);
    data->solver.options().set_reordering_method(strumpack::ReorderingStrategy::METIS);
    data->solver.options().set_compression(strumpack::CompressionType::NONE);
    data->solver.options().set_matching(strumpack::MatchingJob::NONE);
    data->solver.options().disable_replace_tiny_pivots();
    data->solver.options().disable_gpu();
}



template<typename T, typename I>
void solver_csx_strumpack<T,I>::internal_factorize_numeric()
{
    if(data->need_complete) {
        data->op_A_complete.perform_values();
    }

    if(data->called_factorize_numeric) {
        // docs: "The numerical factorization will automatically be redone."
        data->solver.update_matrix_values(data->A_to_use->nrows, data->A_to_use->ptrs, data->A_to_use->idxs, data->A_to_use->vals, is_structurally_symmetric(data->A_to_use->prop.symm));
    }
    else {
        // for some reason I can call reorder only when I have numerical values
        // set_csr_matrix makes and internal copy of the matrix, so maybe thats why
        data->solver.set_csr_matrix(data->A_to_use->nrows, data->A_to_use->ptrs, data->A_to_use->idxs, data->A_to_use->vals, is_structurally_symmetric(data->A_to_use->prop.symm));
        data->solver.reorder();
        data->solver.factor();
    }

    data->called_factorize_numeric = true;
}



template<typename T, typename I>
void solver_csx_strumpack<T,I>::internal_get_permutation(PermutationView_new<I> & perm)
{
    eslog::error("no support for permutation extraction\n");
}



template<typename T, typename I>
void solver_csx_strumpack<T,I>::internal_get_factor_L(MatrixCsxView_new<T,I> & L, bool pattern, bool values)
{
    eslog::error("no support for factor extraction\n");
}



template<typename T, typename I>
void solver_csx_strumpack<T,I>::internal_get_factor_U(MatrixCsxView_new<T,I> & U, bool pattern, bool values)
{
    eslog::error("no support for factor extraction\n");
}



template<typename T, typename I>
void solver_csx_strumpack<T,I>::internal_solve(VectorDenseView_new<T> & rhs, VectorDenseView_new<T> & sol)
{
    MatrixDenseView_new<T> rhs_mat;
    rhs_mat.set_view(rhs.size, 1, rhs.size, 'C', rhs.vals, rhs.ator);

    MatrixDenseView_new<T> sol_mat;
    sol_mat.set_view(sol.size, 1, sol.size, 'C', sol.vals, sol.ator);

    this->internal_solve(rhs_mat, sol_mat);
}



template<typename T, typename I>
void solver_csx_strumpack<T,I>::internal_solve(MatrixDenseView_new<T> & rhs, MatrixDenseView_new<T> & sol)
{
    if(rhs.order != 'C' || sol.order != 'C') eslog::error("only support colmajor rhs/sol for now\n");

    data->solver.solve(rhs.ncols, rhs.vals, rhs.ld, sol.vals, sol.ld);
}



template<typename T, typename I>
void solver_csx_strumpack<T,I>::internal_solve(MatrixCsxView_new<T,I> & rhs, MatrixDenseView_new<T> & sol)
{
    if(sol.order != 'C') eslog::error("only support colmajor sol for now\n");

    MatrixDenseData_new<T> rhs_dn;
    rhs_dn.set(rhs.nrows, rhs.ncols, 'C', AllocatorCPU_new::get_singleton());
    rhs_dn.alloc();

    convert_csx_dny<T,I>::do_all(&rhs, &rhs_dn);

    this->internal_solve(rhs_dn, sol);
}



#define INSTANTIATE_T_I(T,I) \
template class solver_csx_strumpack<T,I>;

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
