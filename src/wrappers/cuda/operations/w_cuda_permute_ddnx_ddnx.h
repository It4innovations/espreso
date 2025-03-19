
#ifndef SRC_WRAPPERS_CUDA_OPERATIONSW_CUDA_PERMUTE_DDNX_DDNX_H
#define SRC_WRAPPERS_CUDA_OPERATIONSW_CUDA_PERMUTE_DDNX_DDNX_H

#include "gpu/operations/permute_ddnx_ddnx.h"



namespace espreso {
namespace gpu {
namespace operations {



#ifdef HAVE_CUDA



template<typename T, typename I>
class w_cuda_permute_ddnx_ddnx : public permute_ddnx_ddnx<T,I>
{
public:
    w_cuda_permute_ddnx_ddnx() = default;
    virtual ~w_cuda_permute_ddnx_ddnx() = default;
protected:
    virtual void internal_setup() override;
    virtual void internal_perform(void * ws_tmp) override;
};



#else



template<typename T, typename I>
class w_cuda_permute_ddnx_ddnx : public submatrix_dcsx_dcsx<T,I>
{
public:
    w_cuda_permute_ddnx_ddnx() { eslog::error("cuda wrapper is not available\n"); }
    virtual ~w_cuda_permute_ddnx_ddnx() = default;
};



#endif



}
}
}

#endif /* SRC_WRAPPERS_CUDA_OPERATIONSW_CUDA_PERMUTE_DDNX_DDNX_H */
