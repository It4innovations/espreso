
#ifndef SRC_WRAPPERS_ONEAPI_OPERATIONS_W_ONEBLAS_HERK_DDNX_DDNY_H
#define SRC_WRAPPERS_ONEAPI_OPERATIONS_W_ONEBLAS_HERK_DDNX_DDNY_H

#include "gpu/operations/herk_ddnx_ddny.h"



namespace espreso {
namespace gpu {
namespace operations {



#ifdef HAVE_ONEAPI



template<typename T>
class w_oneblas_herk_ddnx_ddny : public herk_ddnx_ddny<T>
{
public:
    w_oneblas_herk_ddnx_ddny();
    virtual ~w_oneblas_herk_ddnx_ddny();
protected:
    using herk_ddnx_ddny<T>::q;
    using herk_ddnx_ddny<T>::handle_dnblas;
    using herk_ddnx_ddny<T>::A;
    using herk_ddnx_ddny<T>::C;
    using herk_ddnx_ddny<T>::alpha;
    using herk_ddnx_ddny<T>::beta;
    using herk_ddnx_ddny<T>::mode;
    using herk_ddnx_ddny<T>::wss_tmp_perform;
protected:
    void internal_setup() override;
    void internal_perform(void * /*ws_tmp*/) override;
};



#else



template<typename T>
class w_oneblas_herk_ddnx_ddny : public herk_ddnx_ddny<T>
{
public:
    w_oneblas_herk_ddnx_ddny() { eslog::error("oneapi wrapper is not available\n"); }
    virtual ~w_oneblas_herk_ddnx_ddny() = default;
};



#endif



}
}
}



#endif /* SRC_WRAPPERS_ONEAPI_OPERATIONS_W_ONEBLAS_HERK_DDNX_DDNY_H */
