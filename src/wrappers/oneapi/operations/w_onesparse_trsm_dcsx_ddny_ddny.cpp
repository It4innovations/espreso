
#ifdef HAVE_ONEAPI

#include "wrappers/oneapi/operations/w_onesparse_trsm_dcsx_ddny_ddny.h"

#include "wrappers/oneapi/common_onesparse.h"
#include "wrappers/oneapi/common_oneblas.h"
#include "wrappers/oneapi/common_oneapi_mgm.h"



namespace espreso {
namespace gpu {
namespace operations {



template<typename T, typename I>
w_onesparse_trsm_dcsx_ddny_ddny<T,I>::w_onesparse_trsm_dcsx_ddny_ddny() = default;



template<typename T, typename I>
w_onesparse_trsm_dcsx_ddny_ddny<T,I>::~w_onesparse_trsm_dcsx_ddny_ddny() = default;



template<typename T, typename I>
char w_onesparse_trsm_dcsx_ddny_ddny<T,I>::internal_get_native_place()
{
    return 'O';
}



template<typename T, typename I>
void w_onesparse_trsm_dcsx_ddny_ddny<T,I>::internal_setup()
{
    if(A->order != 'R') eslog::error("A must be rowmajor (onemkl limitation)\n");

    if(place == 'I') {
        ator_ws_tmp_overlap = std::make_unique<AllocatorSinglePointer_new>(AllocatorGPU_new::get_singleton());
        B_tmp.set(B->nrows, B->ncols, B->order, ator_ws_tmp_overlap.get());
    }

    wss_internal = 0; // TODO: check
    wss_persistent = 0;
    wss_tmp_preprocess = 0;
    wss_tmp_perform = ((place == 'I') ? B_tmp.get_memory_impact() : 0);
}



template<typename T, typename I>
void w_onesparse_trsm_dcsx_ddny_ddny<T,I>::internal_preprocess(void * /*ws_tmp*/)
{
    // no-op
}



template<typename T, typename I>
void w_onesparse_trsm_dcsx_ddny_ddny<T,I>::internal_perform(void * ws_tmp)
{
    // my contract allows the user to change sparse matrix values in between calls
    //   this function should then "notice" changes in the values
    // onemkl::sparse contract forbids changing sparse matrix values while they are in the handle
    // so I have to create and destroy the sparse matrix handle with every call
    // furthermore, regarding optimize_trsm:
    // according to some benchmarks, it is not worth it, no change in execution time
    // I need new values in the sparse matrix with every call
    // no update_matrix function, so optimize before every call
    // so it makes no sense to optimize at all

    MatrixDenseView_new<T> * B_to_use = B;

    if(place == 'I') {
        ator_ws_tmp_overlap->set(ws_tmp, wss_tmp_perform);
        B_tmp.alloc();
        B_to_use = &B_tmp;
        oneblas::row_major::omatcopy(q->q, onemkl::transpose::nontrans, B->get_size_primary(), B->get_size_secdary(), T{1}, B->vals, B->ld, B_tmp.vals, B_tmp.ld);
    }

    onesparse::matrix_handle_t handle_A;
    onesparse::init_matrix_handle(&handle_A);
    onesparse::set_csr_data(q->q, handle_A, A->get_size_primary(), A->get_size_secdary(), oneapi::mkl::index_base::zero, A->ptrs, A->idxs, A->vals);

    auto layout = ((X->order == 'R') ? onemkl::layout::row_major : onemkl::layout::col_major);
    auto uplo = ((A->prop.uplo == 'L') ? onemkl::uplo::lower : onemkl::uplo::upper);
    auto diag = ((A->prop.diag == 'U') ? onemkl::diag::unit : onemkl::diag::nonunit);

    onesparse::trsm(q->q, layout, onemkl::transpose::nontrans, onemkl::transpose::nontrans, uplo, diag, T{1}, handle_A, B_to_use->vals, X->ncols, B_to_use->ld, X->vals, X->ld);

    onesparse::release_matrix_handle(q->q, &handle_A);

    if(place == 'I') {
        B_tmp.free();
        ator_ws_tmp_overlap->unset();
    }
}



#define INSTANTIATE_T_I(T,I) \
template class w_onesparse_trsm_dcsx_ddny_ddny<T,I>;

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

#endif
