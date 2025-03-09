
#ifdef HAVE_CUDA
#ifdef USE_CUSPARSE_LEGACY

#include "wrappers/cuda/operations/trsm_dcsx_ddny_ddny.h"

#include <cusparse.h>

#include "wrappers/cuda/common_cusparse.h"



namespace espreso {
namespace gpu {
namespace operations {



template<typename T, typename I>
char w_cusparse_trsm_dcsx_ddny_ddny<T,I>::internal_get_native_place()
{

}



template<typename T, typename I>
void w_cusparse_trsm_dcsx_ddny_ddny<T,I>::internal_set_matrix_A()
{

}



template<typename T, typename I>
void w_cusparse_trsm_dcsx_ddny_ddny<T,I>::internal_set_matrix_X()
{

}



template<typename T, typename I>
void w_cusparse_trsm_dcsx_ddny_ddny<T,I>::internal_set_matrix_B()
{

}



template<typename T, typename I>
void w_cusparse_trsm_dcsx_ddny_ddny<T,I>::internal_setup()
{

}



template<typename T, typename I>
void w_cusparse_trsm_dcsx_ddny_ddny<T,I>::internal_preprocess(void * ws_tmp)
{

}



template<typename T, typename I>
void w_cusparse_trsm_dcsx_ddny_ddny<T,I>::internal_update()
{

}



template<typename T, typename I>
void w_cusparse_trsm_dcsx_ddny_ddny<T,I>::internal_perform(void * ws_tmp)
{

}



#define INSTANTIATE_T_I(T,I) \
template class w_cusparse_trsm_dcsx_ddny_ddny<T,I>;

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
#endif
