
#ifdef HAVE_ONEAPI

#include "wrappers/oneapi/operations/w_oneblas_gemm_ddnx_ddny_ddnz.h"

#include "wrappers/oneapi/common_oneblas.h"
#include "wrappers/oneapi/common_oneapi_mgm.h"



namespace espreso {
namespace gpu {
namespace operations {



template<typename T>
w_oneblas_gemm_ddnx_ddny_ddnz<T>::w_oneblas_gemm_ddnx_ddny_ddnz() = default;



template<typename T>
w_oneblas_gemm_ddnx_ddny_ddnz<T>::~w_oneblas_gemm_ddnx_ddny_ddnz() = default;



template<typename T>
void w_oneblas_gemm_ddnx_ddny_ddnz<T>::internal_setup()
{
    wss_tmp_perform = 0;
}



template<typename T>
void w_oneblas_gemm_ddnx_ddny_ddnz<T>::internal_perform(void * /*ws_tmp*/)
{
    onemkl::transpose op_A = ((A->order == C->order) ? onemkl::transpose::nontrans : onemkl::transpose::trans);
    onemkl::transpose op_B = ((B->order == C->order) ? onemkl::transpose::nontrans : onemkl::transpose::trans);

    if(C->order == 'R') oneblas::row_major   ::gemm(q->q, op_A, op_B, C->nrows, C->ncols, A->ncols, alpha, A->vals, A->ld, B->vals, B->ld, beta, C->vals, C->ld);
    if(C->order == 'C') oneblas::column_major::gemm(q->q, op_A, op_B, C->nrows, C->ncols, A->ncols, alpha, A->vals, A->ld, B->vals, B->ld, beta, C->vals, C->ld);
}



#define INSTANTIATE_T(T) \
template class w_oneblas_gemm_ddnx_ddny_ddnz<T>;

    #define INSTANTIATE \
    /* INSTANTIATE_T(float) */ \
    INSTANTIATE_T(double) \
    /* INSTANTIATE_T(std::complex<float>) */ \
    INSTANTIATE_T(std::complex<double>)

        INSTANTIATE

    #undef INSTANTIATE
#undef INSTANTIATE_T



}
}
}

#endif
