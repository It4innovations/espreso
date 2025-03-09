
#ifndef SRC_WRAPPERS_CUDA_OPERATIONS_HANDLE_CUBLAS_new_H
#define SRC_WRAPPERS_CUDA_OPERATIONS_HANDLE_CUBLAS_new_H

#include "gpu/handle_spblas_new.h"



namespace espreso {
namespace gpu {



#ifdef HAVE_CUDA



struct handle_cusparse_new_internal;

class handle_cusparse_new : public handle_spblas_new
{
public:
    handle_cusparse_new();
    virtual ~handle_cusparse_new();
public:
    std::unique_ptr<handle_cusparse_new_internal> internal;
};



#else



class handle_cusparse_new : public handle_spblas_new
{
    handle_cusparse_new() { eslog::error("cuda wrapper is not available\n"); }
    virtual ~handle_cusparse_new() = default;
};



#endif



}
}

#endif /* SRC_WRAPPERS_CUDA_OPERATIONS_HANDLE_CUBLAS_new_H */
