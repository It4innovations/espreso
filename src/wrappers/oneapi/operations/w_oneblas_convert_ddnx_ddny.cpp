
#ifdef HAVE_ONEAPI

#include "wrappers/oneapi/operations/w_oneblas_convert_ddnx_ddny.h"

#include "wrappers/oneapi/common_oneapi_mgm.h"
#include "wrappers/oneapi/common_oneblas.h"



namespace espreso {
namespace gpu {
namespace operations {



template<typename T>
void w_oneblas_convert_ddnx_ddny<T>::internal_setup()
{
    wss_tmp_perform = 0;
}



template<typename T>
void w_oneblas_convert_ddnx_ddny<T>::internal_perform(void * ws_tmp)
{
    onemkl::transpose op = ((M_src->order == M_dst->order) ? onemkl::transpose::nontrans : onemkl::transpose::trans);

    oneblas::row_major::omatcopy(q->q, op, M_src->get_size_primary(), M_src->get_size_secdary(), T{1}, M_src->vals, M_src->ld, M_dst->vals, M_dst->ld);
}



#define INSTANTIATE_T(T) \
template class w_oneblas_convert_ddnx_ddny<T>;

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
