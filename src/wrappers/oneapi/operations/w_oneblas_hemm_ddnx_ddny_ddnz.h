
#ifndef SRC_WRAPPERS_ONEAPI_OPERATIONS_W_ONEBLAS_HEMM_DDNX_DDNY_DDNZ_H
#define SRC_WRAPPERS_ONEAPI_OPERATIONS_W_ONEBLAS_HEMM_DDNX_DDNY_DDNZ_H

#include "gpu/operations/hemm_ddnx_ddny_ddnz.h"



namespace espreso {
namespace gpu {
namespace operations {



#ifdef HAVE_ONEAPI



template<typename T>
class w_oneblas_hemm_ddnx_ddny_ddnz : public hemm_ddnx_ddny_ddnz<T>
{
public:
    w_oneblas_hemm_ddnx_ddny_ddnz();
    virtual ~w_oneblas_hemm_ddnx_ddny_ddnz();
protected:
    using hemm_ddnx_ddny_ddnz<T>::q;
    using hemm_ddnx_ddny_ddnz<T>::handle_dnblas;
    using hemm_ddnx_ddny_ddnz<T>::A;
    using hemm_ddnx_ddny_ddnz<T>::B;
    using hemm_ddnx_ddny_ddnz<T>::C;
    using hemm_ddnx_ddny_ddnz<T>::wss_tmp_perform;
    using hemm_ddnx_ddny_ddnz<T>::alpha;
    using hemm_ddnx_ddny_ddnz<T>::beta;
protected:
    void internal_setup() override;
    void internal_perform(void * ws_tmp) override;
};



#else



template<typename T>
class w_oneblas_hemm_ddnx_ddny_ddnz : public hemm_ddnx_ddny_ddnz<T>
{
public:
    w_oneblas_hemm_ddnx_ddny_ddnz() { eslog::error("oneapi wrapper is not available\n"); }
    virtual ~w_oneblas_hemm_ddnx_ddny_ddnz() = default;
};



#endif



}
}
}

#endif /* SRC_WRAPPERS_ONEAPI_OPERATIONS_W_ONEBLAS_HEMM_DDNX_DDNY_DDNZ_H */
