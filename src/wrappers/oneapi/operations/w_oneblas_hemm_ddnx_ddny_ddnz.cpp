
#ifdef HAVE_ONEAPI

#include "wrappers/oneapi/operations/w_oneblas_hemm_ddnx_ddny_ddnz.h"

#include "wrappers/oneapi/common_oneblas.h"
#include "wrappers/oneapi/common_oneapi_mgm.h"



namespace espreso {
namespace gpu {
namespace operations {



template<typename T>
w_oneblas_hemm_ddnx_ddny_ddnz<T>::w_oneblas_hemm_ddnx_ddny_ddnz() = default;



template<typename T>
w_oneblas_hemm_ddnx_ddny_ddnz<T>::~w_oneblas_hemm_ddnx_ddny_ddnz() = default;



template<typename T>
void w_oneblas_hemm_ddnx_ddny_ddnz<T>::internal_setup()
{
    wss_tmp_perform = 0;

    if(B->order != C->order) eslog::error("order of matrices B and C must match\n");
    if(utils::is_complex<T>() && A->order != C->order) eslog::error("for complex, order of all matrices must match\n");
}



template<typename T>
void w_oneblas_hemm_ddnx_ddny_ddnz<T>::internal_perform(void * /*ws_tmp*/)
{
    onemkl::uplo uplo = (((C->order == A->order) == (A->prop.uplo == 'L')) ? onemkl::uplo::lower : onemkl::uplo::upper);

    if constexpr(utils::is_complex<T>()) {
        if(C->order == 'R') oneblas::row_major   ::hemm(q->q, onemkl::side::left, uplo, C->nrows, C->ncols, alpha, A->vals, A->ld, B->vals, B->ld, beta, C->vals, C->ld);
        if(C->order == 'C') oneblas::column_major::hemm(q->q, onemkl::side::left, uplo, C->nrows, C->ncols, alpha, A->vals, A->ld, B->vals, B->ld, beta, C->vals, C->ld);
    }
    else {
        if(C->order == 'R') oneblas::row_major   ::symm(q->q, onemkl::side::left, uplo, C->nrows, C->ncols, alpha, A->vals, A->ld, B->vals, B->ld, beta, C->vals, C->ld);
        if(C->order == 'C') oneblas::column_major::symm(q->q, onemkl::side::left, uplo, C->nrows, C->ncols, alpha, A->vals, A->ld, B->vals, B->ld, beta, C->vals, C->ld);
    }
}



#define INSTANTIATE_T(T) \
template class w_oneblas_hemm_ddnx_ddny_ddnz<T>;

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
