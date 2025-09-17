
#ifdef HAVE_MUMPS

#include "wrappers/mumps/operations/schur_csx_dny.mumps.h"

#include "wrappers/mumps/mumps_common.h"
#include "esinfo/mpiinfo.h"
#include "math/primitives_new/allocator_new.h"
#include "math/operations/copy_dnx.h"
#include "math/operations/transpose_dnx_dnx.h"
#include "math/operations/convert_dnx.h"
#include "math/operations/convert_dnx_dny.h"
#include "math/operations/lincomb_matrix_dnx.h"
#include "basis/utilities/stacktimer.h"



namespace espreso {
namespace math {
namespace operations {



template<typename T, typename I>
struct schur_csx_dny_mumps_data
{
    my_mumps_handle_t<T> handle;
    MatrixCsxData_new<T,I> A_whole;
    VectorDenseData_new<I> A_ijv_rowidxs;
    VectorDenseData_new<I> A_ijv_colidxs;
    VectorDenseData_new<I> schur_indices;
    bool is_symmetric = false;
};



template<typename T, typename I>
schur_csx_dny_mumps<T,I>::schur_csx_dny_mumps()
{
    data = std::make_unique<schur_csx_dny_mumps_data<T,I>>();
}



template<typename T, typename I>
schur_csx_dny_mumps<T,I>::~schur_csx_dny_mumps()
{
    data->handle.job = -2; // finalize/terminate
    call_mumps<T>(data->handle);
}



template<typename T, typename I>
void schur_csx_dny_mumps<T,I>::internal_preprocess()
{
    if(called_set_matrix == '4') {
        // todo: I could directly convert the 4 input matrices to the one-based IJV matrix
        this->helper_concat(data->A_whole, 'P');
    }

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

    data->schur_indices.set(size_sc, AllocatorCPU_new::get_singleton());
    data->schur_indices.alloc();
    for(size_t i = 0; i < size_sc; i++) {
        data->schur_indices.vals[i] = size_A11 + i + 1;
    }
    data->handle.size_schur = size_sc;
    data->handle.listvar_schur = data->schur_indices.vals;
    data->handle.icntl[19-1] = 2; // compute schur complement, store only lower triangle (col-major) if symmetric
    data->handle.nprow = 1;
    data->handle.npcol = 1;
    data->handle.mblock = 100;
    data->handle.nblock = 100;

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
void schur_csx_dny_mumps<T,I>::internal_perform_1()
{
    if(called_set_matrix == '4') {
        this->helper_concat(data->A_whole, 'F');
    }

    MatrixDenseData_new<T> sc_tmp; // schur will be in colmajor, if symmetric then in lower
    sc_tmp.set(size_sc, size_sc, 'C', AllocatorCPU_new::get_singleton());
    sc_tmp.alloc();
    sc_tmp.prop.symm = sc->prop.symm;
    if(data->is_symmetric) sc_tmp.prop.uplo = 'L';

    data->handle.a = reinterpret_cast<cpp_type_to_mumps_type_t<T>*>(A->vals);
    data->handle.schur_lld = sc_tmp.ld;
    data->handle.schur = reinterpret_cast<cpp_type_to_mumps_type_t<T>*>(sc_tmp.vals);

    // factorization phase, including schur computation
    data->handle.job = 2;
    call_mumps<T>(data->handle);

    convert_dnx<T>::do_all(&sc_tmp, sc);

    lincomb_matrix_dnx<T>::do_all(sc, alpha, sc, 0, nullptr);
}



template<typename T, typename I>
void schur_csx_dny_mumps<T,I>::internal_perform_2()
{
}



template<typename T, typename I>
void schur_csx_dny_mumps<T,I>::internal_solve_A11(VectorDenseView_new<T> & rhs, VectorDenseView_new<T> & sol)
{
    // mumps needs vector of size of the whole matrix
    // and only considers the elements outside of the schur
    VectorDenseData_new<T> x;
    x.set(size_matrix, AllocatorCPU_new::get_singleton());
    x.alloc();

    std::copy_n(rhs.vals, size_A11, x.vals);

    data->handle.icntl[20-1] = 0; // RHS matrix is dense, and centralized
    data->handle.icntl[21-1] = 0; // solution matrix is centralized
    data->handle.icntl[26-1] = 0; // solve only the internal problem A_11 * x = b;

    // set rhs/sol vector
    data->handle.rhs = reinterpret_cast<cpp_type_to_mumps_type_t<T>*>(x.vals);
    data->handle.nrhs = 1;

    // solve phase
    data->handle.job = 3;
    call_mumps<T>(data->handle);

    std::copy_n(x.vals, size_A11, sol.vals);
}



template<typename T, typename I>
void schur_csx_dny_mumps<T,I>::internal_solve_A11(MatrixDenseView_new<T> & rhs, MatrixDenseView_new<T> & sol)
{
    // not tested

    // mumps needs a rhs/sol matrix as tall as the whole system matrix
    // and only considers the rows outside of the schur
    MatrixDenseData_new<T> tmp;
    tmp.set(size_matrix, rhs.ncols, 'C', AllocatorCPU_new::get_singleton());
    tmp.alloc();

    MatrixDenseView_new<T> tmp_my_sub = tmp.get_submatrix_view(0, size_A11, 0, tmp.ncols);

    math::operations::convert_dnx_dny<T>::do_all(&rhs, &tmp_my_sub, false);

    data->handle.icntl[20-1] = 0; // RHS matrix is dense, and centralized
    data->handle.icntl[21-1] = 0; // solution matrix is centralized
    data->handle.icntl[26-1] = 0; // solve only the internal problem A_11 * x = b;

    // set rhs/sol matrix (col-major)
    data->handle.rhs = reinterpret_cast<cpp_type_to_mumps_type_t<T>*>(tmp.vals);
    data->handle.nrhs = tmp.ncols;
    data->handle.lrhs = tmp.ld;

    // solve phase
    data->handle.job = 3;
    call_mumps<T>(data->handle);

    math::operations::convert_dnx_dny<T>::do_all(&tmp_my_sub, &sol, false);
}



#define INSTANTIATE_T_I(T,I) \
template class schur_csx_dny_mumps<T,I>;

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
