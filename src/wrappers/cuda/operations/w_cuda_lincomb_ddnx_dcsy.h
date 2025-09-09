
#ifndef SRC_WRAPPERS_CUDA_OPERATIONS_W_CUDA_LINCOMB_DDNX_DCSY_H
#define SRC_WRAPPERS_CUDA_OPERATIONS_W_CUDA_LINCOMB_DDNX_DCSY_H

#include "gpu/operations/lincomb_ddnx_dcsy.h"



namespace espreso {
namespace gpu {
namespace operations {



#ifdef HAVE_CUDA



template<typename T, typename I>
class w_cuda_lincomb_ddnx_dcsy : public lincomb_ddnx_dcsy<T,I>
{
public:
    w_cuda_lincomb_ddnx_dcsy() = default;
    virtual ~w_cuda_lincomb_ddnx_dcsy() = default;
protected:
    using lincomb_ddnx_dcsy<T,I>::q;
    using lincomb_ddnx_dcsy<T,I>::X;
    using lincomb_ddnx_dcsy<T,I>::A;
    using lincomb_ddnx_dcsy<T,I>::B;
    using lincomb_ddnx_dcsy<T,I>::alpha;
    using lincomb_ddnx_dcsy<T,I>::beta;
    using lincomb_ddnx_dcsy<T,I>::wss_tmp_perform;
protected:
    virtual void internal_setup() override;
    virtual void internal_perform(void * ws_tmp) override;
};



#else



template<typename T, typename I>
class w_cuda_lincomb_ddnx_dcsy : public lincomb_ddnx_dcsy<T,I>
{
public:
    w_cuda_lincomb_ddnx_dcsy() { eslog::error("cuda wrapper is not available\n"); }
    virtual ~w_cuda_lincomb_ddnx_dcsy() = default;
};



#endif



}
}
}

#endif /* SRC_WRAPPERS_CUDA_OPERATIONS_W_CUDA_LINCOMB_DDNX_DCSY_H */
