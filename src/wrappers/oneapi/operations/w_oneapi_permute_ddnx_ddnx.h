
#ifndef SRC_WRAPPERS_ONEAPI_OPERATIONS_W_ONEAPI_PERMUTE_DDNX_DDNX_H
#define SRC_WRAPPERS_ONEAPI_OPERATIONS_W_ONEAPI_PERMUTE_DDNX_DDNX_H

#include "gpu/operations/permute_ddnx_ddnx.h"



namespace espreso {
namespace gpu {
namespace operations {



#ifdef HAVE_ONEAPI



template<typename T, typename I>
class w_oneapi_permute_ddnx_ddnx : public permute_ddnx_ddnx<T,I>
{
public:
    w_oneapi_permute_ddnx_ddnx() = default;
    virtual ~w_oneapi_permute_ddnx_ddnx() = default;
protected:
    using permute_ddnx_ddnx<T,I>::q;
    using permute_ddnx_ddnx<T,I>::M_src;
    using permute_ddnx_ddnx<T,I>::M_dst;
    using permute_ddnx_ddnx<T,I>::wss_tmp_perform;
    using permute_ddnx_ddnx<T,I>::perm_primary;
    using permute_ddnx_ddnx<T,I>::perm_secdary;
protected:
    void internal_setup() override;
    void internal_perform(void * ws_tmp) override;
};



#else



template<typename T, typename I>
class w_oneapi_permute_ddnx_ddnx : public permute_ddnx_ddnx<T,I>
{
public:
    w_oneapi_permute_ddnx_ddnx() { eslog::error("oneapi wrapper is not available\n"); }
    virtual ~w_oneapi_permute_ddnx_ddnx() = default;
};



#endif



}
}
}

#endif /* SRC_WRAPPERS_ONEAPI_OPERATIONS_W_ONEAPI_PERMUTE_DDNX_DDNX_H */
