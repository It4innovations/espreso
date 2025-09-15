
#ifdef HAVE_ONEAPI

#include "wrappers/oneapi/operations/w_oneblas_herk_ddnx_ddny.h"

#include "wrappers/oneapi/common_oneblas.h"
#include "wrappers/oneapi/common_oneapi_mgm.h"
#include "basis/utilities/utils.h"



namespace espreso {
namespace gpu {
namespace operations {



template<typename T>
w_oneblas_herk_ddnx_ddny<T>::w_oneblas_herk_ddnx_ddny() = default;



template<typename T>
w_oneblas_herk_ddnx_ddny<T>::~w_oneblas_herk_ddnx_ddny() = default;



template<typename T>
void w_oneblas_herk_ddnx_ddny<T>::internal_setup()
{
    if(utils::is_complex<T>()) eslog::error("complex matrices not supported yet\n");

    wss_tmp_perform = 0;
}



template<typename T>
void w_oneblas_herk_ddnx_ddny<T>::internal_perform(void * /*ws_tmp*/)
{
    onemkl::uplo uplo = (((C->prop.uplo == 'L') == (C->order == A->order)) ? onemkl::uplo::lower : onemkl::uplo::upper);
    onemkl::transpose op = ((mode == math::blas::herk_mode::AAh) ? onemkl::transpose::nontrans : onemkl::transpose::trans);
    size_t n = ((mode == math::blas::herk_mode::AAh) ? A->nrows : A->ncols);
    size_t k = ((mode == math::blas::herk_mode::AAh) ? A->ncols : A->nrows);

    if constexpr(utils::is_complex<T>()) {
        eslog::error("complex matrices not supported yet\n");
    }
    else {
        if(A->order == 'R') oneblas::row_major   ::syrk(q->q, uplo, op, n, k, alpha, A->vals, A->ld, beta, C->vals, C->ld);
        if(A->order == 'C') oneblas::column_major::syrk(q->q, uplo, op, n, k, alpha, A->vals, A->ld, beta, C->vals, C->ld);
    }
}



#define INSTANTIATE_T(T) \
template class w_oneblas_herk_ddnx_ddny<T>;

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
