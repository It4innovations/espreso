
#ifdef HAVE_PASTIX

#include "wrappers/pastix/operations/schur_csx_dny.pastix.h"

#include "esinfo/ecfinfo.h"
#include "wrappers/pastix/pastix_common.h"
#include "math/primitives_new/allocator_new.h"
#include "math/operations/convert_dnx.h"
#include "math/operations/lincomb_matrix_dnx.h"



namespace espreso {
namespace math {
namespace operations {



template<typename T, typename I>
struct schur_csx_dny_pastix_data
{
    MatrixCsxData_new<T,I> A_whole;
    VectorDenseData_new<pastix_int_t> schur_unknowns;
    pastix_int_t iparm[IPARM_SIZE];
    double dparm[DPARM_SIZE];
    pastix_data_t * pastix_data = nullptr;
    MatrixCsxView_new<T,I> A_to_use;
    spmatrix_t pastix_A;
    struct config {
        bool use_gpu = false;
    } cfg;
    bool cholesky;
};



template<typename T, typename I>
static void setup_config(typename schur_csx_dny_pastix_data<T,I>::config & cfg)
{
    using ecf_config = SchurCsxDnyPastixConfig;
    const ecf_config & ecf = info::ecf->operations.schur_csx_dny_pastix;

    switch(ecf.use_gpu) {
        case ecf_config::AUTOBOOL::AUTO:  cfg.use_gpu = false; break;
        case ecf_config::AUTOBOOL::TRUE:  cfg.use_gpu = true;  break;
        case ecf_config::AUTOBOOL::FALSE: cfg.use_gpu = false; break;
    }
}



template<typename T, typename I>
schur_csx_dny_pastix<T,I>::schur_csx_dny_pastix()
{
    eslog::error("schur_csx_dny_pastix produces wrong results and should not be used. Feel free to comment out this error to test it\n");

    data = std::make_unique<schur_csx_dny_pastix_data<T,I>>();
    setup_config<T,I>(data->cfg);

    check_pastix_instances(data->cfg.use_gpu, true);
}



template<typename T, typename I>
schur_csx_dny_pastix<T,I>::~schur_csx_dny_pastix()
{
    if(data->pastix_data != nullptr) {
        pastixFinalize(&data->pastix_data);
    }

    check_pastix_instances(data->cfg.use_gpu, false);
}



template<typename T, typename I>
void schur_csx_dny_pastix<T,I>::internal_preprocess()
{
    if(called_set_matrix == '4') {
        this->helper_concat(data->A_whole, 'P');
    }

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
    data->iparm[IPARM_SCHUR_FACT_MODE] = PastixFactModeLocal;
    data->iparm[IPARM_SCHUR_SOLV_MODE] = PastixSolvModeLocal;

    pastixInit(&data->pastix_data, 0, data->iparm, data->dparm);

    if(A->order == 'C') {
        data->A_to_use = *A;
    }
    else {
        data->A_to_use = A->get_transposed_reordered_view();
    }
    
    spmInit(&data->pastix_A);
    data->pastix_A.mtxtype = symm_to_pastix<T>(data->A_to_use.prop.symm);
    data->pastix_A.flttype = type_to_pastix<T>();
    data->pastix_A.fmttype = order_to_pastix(data->A_to_use.order);
    data->pastix_A.baseval = 0;
    data->pastix_A.n = data->A_to_use.nrows;
    data->pastix_A.nnz = data->A_to_use.nnz;
    data->pastix_A.dof = 1;
    data->pastix_A.colptr = data->A_to_use.ptrs; // csc
    data->pastix_A.rowptr = data->A_to_use.idxs; // csc
    data->pastix_A.values = data->A_to_use.vals;
    data->pastix_A.replicated = -1; // neither distrubuted nor replicated, just one matrix
    spmUpdateComputedFields(&data->pastix_A);

    data->schur_unknowns.set(size_sc, AllocatorCPU_new::get_singleton());
    data->schur_unknowns.alloc();
    for(size_t i = 0; i < size_sc; i++) {
        data->schur_unknowns.vals[i] = size_A11 + i;
    }
    pastixSetSchurUnknownList(data->pastix_data, size_sc, data->schur_unknowns.vals);

    pastix_task_analyze(data->pastix_data, &data->pastix_A);
}



template<typename T, typename I>
void schur_csx_dny_pastix<T,I>::internal_perform_1()
{
    if(called_set_matrix == '4') {
        this->helper_concat(data->A_whole, 'F');
    }

    if(A->order == 'C') {
        data->A_to_use = *A;
    }
    else {
        data->A_to_use = A->get_transposed_reordered_view();
    }

    data->pastix_A.colptr = data->A_to_use.ptrs; // csc
    data->pastix_A.rowptr = data->A_to_use.idxs; // csc
    data->pastix_A.values = data->A_to_use.vals;

    pastix_task_numfact(data->pastix_data, &data->pastix_A);
}



template<typename T, typename I>
void schur_csx_dny_pastix<T,I>::internal_perform_2()
{
    MatrixDenseData_new<T> sc_tmp;
    sc_tmp.set(size_sc, size_sc, 'C', AllocatorCPU_new::get_singleton());
    sc_tmp.prop.symm = sc->prop.symm;
    sc_tmp.alloc();

    pastixGetSchur(data->pastix_data, sc_tmp.vals, sc_tmp.ld);

    convert_dnx<T>::do_all(&sc_tmp, sc);

    if(alpha != T{1}) {
        lincomb_matrix_dnx<T>::do_all(sc, alpha, sc, 0, nullptr);
    }
}



template<typename T, typename I>
void schur_csx_dny_pastix<T,I>::internal_solve_A11(VectorDenseView_new<T> & rhs, VectorDenseView_new<T> & sol)
{
    if(&rhs != &sol) {
        std::copy_n(rhs.vals, rhs.size, sol.vals);
    }

    pastix_task_solve(data->pastix_data, data->pastix_A.gNexp, 1, sol.vals, size_A11);
}



#define INSTANTIATE_T_I(T,I) \
template class schur_csx_dny_pastix<T,I>;

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
