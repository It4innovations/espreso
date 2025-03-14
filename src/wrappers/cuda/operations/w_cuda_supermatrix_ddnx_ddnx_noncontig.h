
#ifndef SRC_WRAPPERS_CUDA_OPERATIONS_W_CUDA_SUPERMATRIX_DDNX_DDNX_H
#define SRC_WRAPPERS_CUDA_OPERATIONS_W_CUDA_SUPERMATRIX_DDNX_DDNX_H

#include "gpu/operations/supermatrix_ddnx_ddnx_noncontig.h"



namespace espreso {
namespace gpu {
namespace operations {



#ifdef HAVE_CUDA



template<typename T>
class w_cuda_supermatrix_ddnx_ddnx_noncontig : public supermatrix_ddnx_ddnx_noncontig<T>
{
public:
    w_cuda_supermatrix_ddnx_ddnx_noncontig() = default;
    virtual ~w_cuda_supermatrix_ddnx_ddnx_noncontig() = default;
protected:
    virtual void internal_perform() override;
};



#else



template<typename T>
class w_cuda_supermatrix_ddnx_ddnx_noncontig : public supermatrix_ddnx_ddnx_noncontig<T>
{
public:
    w_cuda_supermatrix_ddnx_ddnx_noncontig() { eslog::error("cuda wrapper is not available\n"); }
    virtual ~w_cuda_supermatrix_ddnx_ddnx_noncontig() = default;
};



#endif



}
}
}

#endif /* SRC_WRAPPERS_CUDA_OPERATIONS_W_CUDA_SUPERMATRIX_DDNX_DDNX_H */
