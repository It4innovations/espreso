
#ifdef HAVE_ONEAPI

#include "wrappers/oneapi/operations/w_onesparse_gemm_dcsx_ddny_ddnz.h"

#include "wrappers/oneapi/common_onesparse.h"



namespace espreso {
namespace gpu {
namespace operations {



template<typename T, typename I>
w_onesparse_gemm_dcsx_ddny_ddnz<T,I>::w_onesparse_gemm_dcsx_ddny_ddnz() = default;



template<typename T, typename I>
w_onesparse_gemm_dcsx_ddny_ddnz<T,I>::~w_onesparse_gemm_dcsx_ddny_ddnz() = default;



template<typename T, typename I>
void w_onesparse_gemm_dcsx_ddny_ddnz<T,I>::internal_setup()
{
    // onemkl::sparse currently supports only nontrans for matrix B
    if(B->order != C->order) eslog::error("B and C must have equal order\n");

    wss_internal = 0; // TODO: check
    wss_persistent = 0;
    wss_tmp_preprocess = 0;
    wss_tmp_perform = 0;
}



template<typename T, typename I>
void w_onesparse_gemm_dcsx_ddny_ddnz<T,I>::internal_preprocess(void * /*ws_tmp*/)
{
    // no-op
}



template<typename T, typename I>
void w_onesparse_gemm_dcsx_ddny_ddnz<T,I>::internal_perform(void * /*ws_tmp*/)
{
    // my contract allows the user to change sparse matrix values in between calls
    //   this function should then "notice" changes in the values
    // onemkl::sparse contract forbids changing sparse matrix values while they are in the handle
    // so I have to create and destroy the sparse matrix handle with every call
    // regarding optimize_gemm:
    // I need new values in the sparse matrix with every call
    // no update_matrix function, so optimize before every call
    // so it makes no sense to optimize at all

    onesparse::matrix_handle_t handle_A;

    onesparse::init_matrix_handle(&handle_A);
    onesparse::set_csr_data(q->q, handle_A, A->get_size_primary(), A->get_size_secdary(), oneapi::mkl::index_base::zero, A->ptrs, A->idxs, A->vals);

    auto layout = ((C->order == 'R') ? onemkl::layout::row_major : onemkl::layout::col_major);
    auto op_A = ((A->order == 'R') ? onemkl::transpose::nontrans : onemkl::transpose::trans);
    auto op_B = ((B->order == C->order) ? onemkl::transpose::nontrans : onemkl::transpose::trans);

    onesparse::gemm(q->q, layout, op_A, op_B, alpha, handle_A, B->vals, B->ncols, B->ld, beta, C->vals, C->ld);

    onesparse::release_matrix_handle(q->q, &handle_A);
}



#define INSTANTIATE_T_I(T,I) \
template class w_onesparse_gemm_dcsx_ddny_ddnz<T,I>;

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
