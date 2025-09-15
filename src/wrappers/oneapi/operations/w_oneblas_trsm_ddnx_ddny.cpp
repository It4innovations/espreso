
#ifdef HAVE_ONEAPI

#include "wrappers/oneapi/operations/w_oneblas_trsm_ddnx_ddny.h"

#include "wrappers/oneapi/common_oneblas.h"
#include "wrappers/oneapi/common_oneapi_mgm.h"



namespace espreso {
namespace gpu {
namespace operations {



template<typename T>
w_oneblas_trsm_ddnx_ddny<T>::w_oneblas_trsm_ddnx_ddny() = default;



template<typename T>
w_oneblas_trsm_ddnx_ddny<T>::~w_oneblas_trsm_ddnx_ddny() = default;



template<typename T>
void w_oneblas_trsm_ddnx_ddny<T>::internal_setup()
{
    wss_tmp_perform = 0;
}



template<typename T>
void w_oneblas_trsm_ddnx_ddny<T>::internal_perform(void * ws_tmp)
{
    onemkl::uplo uplo = (((A->order == X->order) == (A->prop.uplo == 'L')) ? onemkl::uplo::lower : onemkl::uplo::upper);
    onemkl::transpose op = ((A->order == X->order) ? onemkl::transpose::nontrans : onemkl::transpose::trans);
    onemkl::diag diag = ((A->prop.diag == 'U') ? onemkl::diag::unit : onemkl::diag::nonunit);

    if(X->order == 'R') oneblas::row_major   ::trsm(q->q, onemkl::side::left, uplo, op, diag, X->nrows, X->ncols, T{1}, A->vals, A->ld, X->vals, X->ld);
    if(X->order == 'C') oneblas::column_major::trsm(q->q, onemkl::side::left, uplo, op, diag, X->nrows, X->ncols, T{1}, A->vals, A->ld, X->vals, X->ld);
}



#define INSTANTIATE_T(T) \
template class w_oneblas_trsm_ddnx_ddny<T>;

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
