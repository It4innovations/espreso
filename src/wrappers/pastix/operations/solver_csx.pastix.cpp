
#ifdef HAVE_PASTIX

#include "wrappers/pastix/operations/solver_csx.pastix.h"

#include "wrappers/pastix/pastix_common.h"
#include "math/primitives_new/allocator_new.h"
#include "math/operations/convert_dnx_dny.h"
#include "math/operations/convert_csx_dny.h"
#include "math/operations/convert_csx_csy.h"
#include "math/operations/copy_dnx.h"



namespace espreso {
namespace math {
namespace operations {



template<typename T, typename I>
struct solver_csx_pastix_data
{
    pastix_int_t iparm[IPARM_SIZE];
    double dparm[DPARM_SIZE];
    pastix_data_t * pastix_data = nullptr;
    spmatrix_t pastix_A;
    struct config {
        bool use_gpu = false;
    } cfg;
    bool cholesky;
};



template<typename T, typename I>
solver_csx_pastix<T,I>::solver_csx_pastix()
{
    data = std::make_unique<solver_csx_pastix_data<T,I>>();

    check_pastix_instances(data->cfg.use_gpu, true);
}



template<typename T, typename I>
solver_csx_pastix<T,I>::~solver_csx_pastix()
{
    check_pastix_instances(data->cfg.use_gpu, false);
}



template<typename T, typename I>
void solver_csx_pastix<T,I>::internal_factorize_symbolic()
{
    data->cholesky = is_hermitian<T>(A->prop.symm);

    pastixInitParam(data->iparm, data->dparm);

    data->iparm[IPARM_THREAD_NBR] = 1;
    data->iparm[IPARM_VERBOSE] = PastixVerboseNot;
    data->iparm[IPARM_ORDERING] = PastixOrderMetis;
    data->iparm[IPARM_FACTORIZATION] = (data->cholesky ? PastixFactLLH : PastixFactLU);
    if(data->cfg.use_gpu) {
        data->iparm[IPARM_GPU_NBR] = 1;
        data->iparm[IPARM_SCHEDULER] = PastixSchedParsec;
    }
    else {
        data->iparm[IPARM_SCHEDULER] = PastixSchedSequential;
    }

    pastixInit(&data->pastix_data, 0, data->iparm, data->dparm);

    spmInit(&data->pastix_A);
    data->pastix_A.mtxtype = symm_to_pastix<T>(A->prop.symm);
    data->pastix_A.flttype = type_to_pastix<T>();
    data->pastix_A.fmttype = order_to_pastix(A->order);
    data->pastix_A.baseval = 0;
    data->pastix_A.n = A->nrows;
    data->pastix_A.nnz = A->nnz;
    data->pastix_A.dof = 1;
    data->pastix_A.rowptr = ((A->order == 'R') ? A->ptrs : A->idxs);
    data->pastix_A.colptr = ((A->order == 'R') ? A->idxs : A->ptrs);
    data->pastix_A.values = A->vals;
    data->pastix_A.replicated = -1; // neither distrubuted nor replicated, just one matrix
    spmUpdateComputedFields(&data->pastix_A);

    pastix_task_analyze(data->pastix_data, &data->pastix_A);
}



template<typename T, typename I>
void solver_csx_pastix<T,I>::internal_factorize_numeric()
{
    data->pastix_A.colptr = A->ptrs;
    data->pastix_A.rowptr = A->idxs;
    data->pastix_A.values = A->vals;

    pastix_task_numfact(data->pastix_data, &data->pastix_A);
}



template<typename T, typename I>
void solver_csx_pastix<T,I>::internal_get_permutation(PermutationView_new<I> & perm)
{
    eslog::error("no support for permutation extraction\n");
}



template<typename T, typename I>
void solver_csx_pastix<T,I>::internal_get_factor_L(MatrixCsxView_new<T,I> & L, bool pattern, bool values)
{
    eslog::error("no support for factor extraction\n");
}



template<typename T, typename I>
void solver_csx_pastix<T,I>::internal_get_factor_U(MatrixCsxView_new<T,I> & U, bool pattern, bool values)
{
    eslog::error("no support for factor extraction\n");
}



template<typename T, typename I>
void solver_csx_pastix<T,I>::internal_solve(VectorDenseView_new<T> & rhs, VectorDenseView_new<T> & sol)
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
void solver_csx_pastix<T,I>::internal_solve(MatrixDenseView_new<T> & rhs, MatrixDenseView_new<T> & sol)
{
    if(sol.order == 'R') {
        MatrixDenseData_new<T> tmp;
        tmp.set(rhs.nrows, rhs.ncols, 'C', AllocatorCPU_new::get_singleton());
        tmp.alloc();
        convert_dnx_dny<T>::do_all(&rhs, &tmp, false);
        this->internal_solve(tmp, tmp);
        convert_dnx_dny<T>::do_all(&tmp, &sol, false);
        return;
    }

    if(&rhs != &sol) {
        copy_dnx<T>::do_all(&rhs, &sol, false);
    }

    pastix_task_solve(data->pastix_data, data->pastix_A.gNexp, sol.ncols, sol.vals, sol.ld);
}



template<typename T, typename I>
void solver_csx_pastix<T,I>::internal_solve(MatrixCsxView_new<T,I> & rhs, MatrixDenseView_new<T> & sol)
{
    if(sol.order == 'R') {
        MatrixDenseData_new<T> tmp;
        tmp.set(rhs.nrows, rhs.ncols, 'C', AllocatorCPU_new::get_singleton());
        tmp.alloc();
        convert_csx_dny<T,I>::do_all(&rhs, &tmp);
        this->internal_solve(tmp, tmp);
        convert_dnx_dny<T>::do_all(&tmp, &sol, false);
    }
    else {
        convert_csx_dny<T,I>::do_all(&rhs, &sol);
        this->internal_solve(sol, sol);
    }
}



#define INSTANTIATE_T_I(T,I) \
template class solver_csx_pastix<T,I>;

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
