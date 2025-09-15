
#ifndef SRC_WRAPPERS_ONEAPI_OPERATIONSW_ONEBLAS_COPY_DDNX_DDNX_H
#define SRC_WRAPPERS_ONEAPI_OPERATIONSW_ONEBLAS_COPY_DDNX_DDNX_H

#include "gpu/operations/convert_ddnx_ddny.h"



namespace espreso {
namespace gpu {
namespace operations {



#ifdef HAVE_ONEAPI



template<typename T>
class w_oneblas_convert_ddnx_ddny : public convert_ddnx_ddny<T>
{
public:
    w_oneblas_convert_ddnx_ddny() = default;
    virtual ~w_oneblas_convert_ddnx_ddny() = default;
protected:
    using convert_ddnx_ddny<T>::q;
    using convert_ddnx_ddny<T>::handle_dnblas;
    using convert_ddnx_ddny<T>::M_src;
    using convert_ddnx_ddny<T>::M_dst;
    using convert_ddnx_ddny<T>::wss_tmp_perform;
protected:
    void internal_setup() override;
    void internal_perform(void * ws_tmp) override;
};



#else



template<typename T>
class w_oneblas_convert_ddnx_ddny : public convert_ddnx_ddny<T>
{
public:
    w_oneblas_convert_ddnx_ddny() { eslog::error("oneapi wrapper is not available\n"); }
    virtual ~w_oneblas_convert_ddnx_ddny() = default;
};



#endif



}
}
}

#endif /* SRC_WRAPPERS_ONEAPI_OPERATIONSW_ONEBLAS_COPY_DDNX_DDNX_H */
