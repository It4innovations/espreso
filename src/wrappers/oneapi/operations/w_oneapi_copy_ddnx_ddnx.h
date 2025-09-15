
#ifndef SRC_WRAPPERS_ONEAPI_OPERATIONS_W_ONEAPI_COPY_DDNX_DDNX_H
#define SRC_WRAPPERS_ONEAPI_OPERATIONS_W_ONEAPI_COPY_DDNX_DDNX_H

#include "gpu/operations/copy_ddnx_ddnx.h"



namespace espreso {
namespace gpu {
namespace operations {



#ifdef HAVE_ONEAPI



template<typename T>
class w_oneapi_copy_ddnx_ddnx : public copy_ddnx_ddnx<T>
{
public:
    w_oneapi_copy_ddnx_ddnx() = default;
    virtual ~w_oneapi_copy_ddnx_ddnx() = default;
protected:
    using copy_ddnx_ddnx<T>::q;
    using copy_ddnx_ddnx<T>::M_src;
    using copy_ddnx_ddnx<T>::M_dst;
    using copy_ddnx_ddnx<T>::uplo;
    using copy_ddnx_ddnx<T>::wss_tmp_perform;
protected:
    void internal_setup() override;
    void internal_perform(void * ws_tmp) override;
};



#else



template<typename T>
class w_oneapi_copy_ddnx_ddnx : public copy_ddnx_ddnx<T>
{
public:
    w_oneapi_copy_ddnx_ddnx() { eslog::error("oneapi wrapper is not available\n"); }
    virtual ~w_oneapi_copy_ddnx_ddnx() = default;
};



#endif



}
}
}

#endif /* SRC_WRAPPERS_ONEAPI_OPERATIONS_W_ONEAPI_COPY_DDNX_DDNX_H */
