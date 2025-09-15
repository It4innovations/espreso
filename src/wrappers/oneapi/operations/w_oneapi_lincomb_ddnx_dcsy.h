
#ifndef SRC_WRAPPERS_ONEAPI_OPERATIONS_W_ONEAPI_LINCOMB_DDNX_DCSY_H
#define SRC_WRAPPERS_ONEAPI_OPERATIONS_W_ONEAPI_LINCOMB_DDNX_DCSY_H

#include "gpu/operations/lincomb_ddnx_dcsy.h"



namespace espreso {
namespace gpu {
namespace operations {



#ifdef HAVE_ONEAPI



template<typename T, typename I>
class w_oneapi_lincomb_ddnx_dcsy : public lincomb_ddnx_dcsy<T,I>
{
public:
    w_oneapi_lincomb_ddnx_dcsy() = default;
    virtual ~w_oneapi_lincomb_ddnx_dcsy() = default;
protected:
    using lincomb_ddnx_dcsy<T,I>::q;
    using lincomb_ddnx_dcsy<T,I>::X;
    using lincomb_ddnx_dcsy<T,I>::A;
    using lincomb_ddnx_dcsy<T,I>::B;
    using lincomb_ddnx_dcsy<T,I>::alpha;
    using lincomb_ddnx_dcsy<T,I>::beta;
    using lincomb_ddnx_dcsy<T,I>::wss_tmp_perform;
protected:
    void internal_setup() override;
    void internal_perform(void * ws_tmp) override;
};



#else



template<typename T, typename I>
class w_oneapi_lincomb_ddnx_dcsy : public lincomb_ddnx_dcsy<T,I>
{
public:
    w_oneapi_lincomb_ddnx_dcsy() { eslog::error("oneapi wrapper is not available\n"); }
    virtual ~w_oneapi_lincomb_ddnx_dcsy() = default;
};



#endif



}
}
}

#endif /* SRC_WRAPPERS_ONEAPI_OPERATIONS_W_ONEAPI_LINCOMB_DDNX_DCSY_H */
