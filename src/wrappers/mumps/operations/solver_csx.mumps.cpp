
#ifdef HAVE_MUMPS

#include "wrappers/mumps/operations/solver_csx.mumps.h"

#include "wrappers/mumps/mumps_common.h"
#include "esinfo/mpiinfo.h"
#include "math/primitives_new/allocator_new.h"
#include "math/operations/copy_dnx.h"
#include "math/operations/convert_dnx_dny.h"
#include "math/operations/convert_csx_csy.h"



namespace espreso {
namespace math {
namespace operations {



template<typename T, typename I>
struct solver_csx_mumps_data
{
    my_mumps_handle_t<T> handle;
    VectorDenseData_new<I> A_ijv_rowidxs;
    VectorDenseData_new<I> A_ijv_colidxs;
    bool is_symmetric = false;
};



template<typename T, typename I>
solver_csx_mumps<T,I>::solver_csx_mumps()
{
    data = std::make_unique<solver_csx_mumps_data<T,I>>();
}



template<typename T, typename I>
solver_csx_mumps<T,I>::~solver_csx_mumps()
{
    data->handle.job = -2; // finalize/terminate
    call_mumps<T>(data->handle);
}



template<typename T, typename I>
void solver_csx_mumps<T,I>::internal_factorize_symbolic()
{
    // mumps handles hermitian matrices as general nonsymmetric. TODO figure it out
    // btw, hermitian matrices must have no zeros on diagonal, so no problem figuring out nnz of the full matrix
    if constexpr(utils::is_complex<T>()) if(is_hermitian<T>(A->prop.symm) && is_uplo(A->prop.uplo)) eslog::error("no support for complex hermitian matrices yet\n");

    // convert csx matrix to one-based IJV. Only the row+col indices are required
    mumps_helper_csx_to_ijv<T,I>(*A, data->A_ijv_rowidxs, data->A_ijv_colidxs);

    data->is_symmetric = is_symmetric<T>(A->prop.symm);

    data->handle.sym = (data->is_symmetric ? 1 : 0); // symmetric or unsymmetric matrix
    data->handle.par = 1; // this process is involved in computations
    data->handle.comm_fortran = info::mpi::comm_c2f(MPI_COMM_SELF);

    data->handle.job = -1; // initialize
    call_mumps<T>(data->handle);

    data->handle.icntl[4-1] = 0; // print only error messages, no statistics
    data->handle.icntl[6-1] = 0; // no permutation based on numerical values on diagonal
    data->handle.icntl[8-1] = 0; // no row/col scaling
    data->handle.icntl[28-1] = 1; // ordering will be computed sequentially
    data->handle.icntl[7-1] = 5; // use metis for fill-reducing ordering
    data->handle.icntl[10-1] = 0; // no iterative refinement
    data->handle.icntl[16-1] = 1; // use only a single OpenMP thread (mumps calls omp_set_num_thread inside with this number)

    // input sparse matrix in IJV/COO format
    data->handle.icntl[5-1] = 0; // matrix is in IJV/COO format
    data->handle.icntl[18-1] = 0; // the whole matrix is present here (centralized, no distribution)
    data->handle.n = A->nrows;
    data->handle.nnz = A->nnz;
    data->handle.irn = data->A_ijv_rowidxs.vals;
    data->handle.jcn = data->A_ijv_colidxs.vals;

    // analysis phase
    data->handle.job = 1;
    call_mumps<T>(data->handle);
}



template<typename T, typename I>
void solver_csx_mumps<T,I>::internal_factorize_numeric()
{
    data->handle.a = reinterpret_cast<cpp_type_to_mumps_type_t<T>*>(A->vals);

    // factorization phase
    data->handle.job = 2;
    call_mumps<T>(data->handle);
}



template<typename T, typename I>
void solver_csx_mumps<T,I>::internal_get_permutation(PermutationView_new<I> & perm)
{
    if constexpr(std::is_same_v<I,int>) {
        std::copy_n(data->handle.sym_perm, perm.size, perm.src_to_dst);
    }
    else {
        auto begin = data->handle.sym_perm;
        auto end = begin + perm.size;
        std::transform(begin, end, perm.src_to_dst, [](int i){ return static_cast<I>(i); });
    }

    PermutationView_new<I>::invert(perm.src_to_dst, perm.dst_to_src, perm.size);
}



template<typename T, typename I>
void solver_csx_mumps<T,I>::internal_get_factor_L(MatrixCsxView_new<T,I> & L, bool pattern, bool values)
{
    eslog::error("no support for factor extraction\n");
}



template<typename T, typename I>
void solver_csx_mumps<T,I>::internal_get_factor_U(MatrixCsxView_new<T,I> & U, bool pattern, bool values)
{
    eslog::error("no support for factor extraction\n");
}



template<typename T, typename I>
void solver_csx_mumps<T,I>::internal_solve(VectorDenseView_new<T> & rhs, VectorDenseView_new<T> & sol)
{
    MatrixDenseView_new<T> rhs_mat;
    rhs_mat.set_view(rhs.size, 1, rhs.size, 'C', rhs.vals, rhs.ator);

    MatrixDenseView_new<T> sol_mat;
    sol_mat.set_view(sol.size, 1, sol.size, 'C', sol.vals, sol.ator);

    if(&rhs == &sol) {
        this->internal_solve(sol_mat, sol_mat);
    }
    else {
        this->internal_solve(rhs_mat, sol_mat);
    }
}



template<typename T, typename I>
void solver_csx_mumps<T,I>::internal_solve(MatrixDenseView_new<T> & rhs, MatrixDenseView_new<T> & sol)
{
    if(&rhs != &sol) {
        copy_dnx<T>::do_all(&rhs, &sol, false);
    }

    if(sol.order != 'C') {
        MatrixDenseData_new<T> tmp;
        tmp.set(sol.ncols, sol.nrows, 'C', AllocatorCPU_new::get_singleton());
        tmp.alloc();
        convert_dnx_dny<T>::do_all(&sol, &tmp, false);
        this->internal_solve(tmp, tmp);
        convert_dnx_dny<T>::do_all(&tmp, &sol, false);
        return;
    }
    
    data->handle.icntl[20-1] = 0; // RHS matrix is dense, and centralized
    data->handle.icntl[21-1] = 0; // solution matrix is centralized

    // set rhs/sol matrix
    data->handle.rhs = reinterpret_cast<cpp_type_to_mumps_type_t<T>*>(sol.vals);
    data->handle.nrhs = sol.ncols;
    data->handle.lrhs = sol.ld; // leading dimension

    // solve phase
    data->handle.job = 3;
    call_mumps<T>(data->handle);
}



template<typename T, typename I>
void solver_csx_mumps<T,I>::internal_solve(MatrixCsxView_new<T,I> & rhs, MatrixDenseView_new<T> & sol)
{
    if(rhs.order != 'C') {
        MatrixCsxData_new<T,I> rhs_2;
        rhs_2.set(rhs.nrows, rhs.ncols, rhs.nnz, 'C', AllocatorCPU_new::get_singleton());
        rhs_2.alloc();
        convert_csx_csy<T,I>::do_all(&rhs, &rhs_2);
        this->internal_solve(rhs_2, sol);
        return;
    }

    if(sol.order != 'C') {
        MatrixDenseData_new<T> sol_2;
        sol_2.set(sol.ncols, sol.nrows, 'C', AllocatorCPU_new::get_singleton());
        sol_2.alloc();
        this->internal_solve(rhs, sol_2);
        convert_dnx_dny<T>::do_all(&sol_2, &sol, false);
        return;
    }

    // mumps requires one-based indexing. Copy the rhs matrix and add 1 to the ptrs and idxs
    VectorDenseData_new<I> rhs_onebased_ptrs;
    VectorDenseData_new<I> rhs_onebased_idxs;
    rhs_onebased_ptrs.set(rhs.get_size_primary() + 1, AllocatorCPU_new::get_singleton());
    rhs_onebased_idxs.set(rhs.nnz, AllocatorCPU_new::get_singleton());
    rhs_onebased_ptrs.alloc();
    rhs_onebased_idxs.alloc();

    auto add_one = [](I i){ return i + 1; };
    utils::transform_n(rhs.ptrs, rhs.get_size_primary() + 1, rhs_onebased_ptrs.vals, add_one);
    utils::transform_n(rhs.idxs, rhs.nnz, rhs_onebased_idxs.vals, add_one);

    data->handle.icntl[20-1] = 3; // RHS matrix is sparse centralized, and the sparsity is utilized
    data->handle.icntl[21-1] = 0; // solution matrix is centralized

    data->handle.nz_rhs = rhs.nnz;
    data->handle.nrhs = rhs.ncols;
    data->handle.irhs_ptr = rhs_onebased_ptrs.vals;
    data->handle.irhs_sparse = rhs_onebased_idxs.vals;
    data->handle.rhs_sparse = reinterpret_cast<cpp_type_to_mumps_type_t<T>*>(rhs.vals);
    data->handle.rhs = reinterpret_cast<cpp_type_to_mumps_type_t<T>*>(sol.vals);
    data->handle.lrhs = sol.ld;

    // solve phase
    data->handle.job = 3;
    call_mumps<T>(data->handle);
}



#define INSTANTIATE_T_I(T,I) \
template class solver_csx_mumps<T,I>;

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
