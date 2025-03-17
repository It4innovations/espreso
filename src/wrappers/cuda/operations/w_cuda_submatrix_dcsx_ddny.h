
#ifndef SRC_WRAPPERS_CUDA_OPERATIONS_W_CUDA_SUBMATRIX_DCSX_DDNY_H
#define SRC_WRAPPERS_CUDA_OPERATIONS_W_CUDA_SUBMATRIX_DCSX_DDNY_H

#include "gpu/operations/submatrix_dcsx_ddny.h"



namespace espreso {
namespace gpu {
namespace operations {



#ifdef HAVE_CUDA



template<typename T, typename I>
class w_cuda_submatrix_dcsx_ddny : public submatrix_dcsx_ddny<T,I>
{
public:
    w_cuda_submatrix_dcsx_ddny() = default;
    virtual ~w_cuda_submatrix_dcsx_ddny() = default;
protected:
    virtual void internal_setup() override;
    virtual void internal_preprocess(void * ws_tmp) override;
    virtual void internal_perform(void * ws_tmp) override;
};



#else



template<typename T, typename I>
class w_cuda_submatrix_dcsx_ddny : public submatrix_dcsx_ddny<T,I>
{
public:
    w_cuda_submatrix_dcsx_ddny() { eslog::error("cuda wrapper is not available\n"); }
    virtual ~w_cuda_submatrix_dcsx_ddny() = default;
};



#endif



}
}
}

#endif /* SRC_WRAPPERS_CUDA_OPERATIONS_W_CUDA_SUBMATRIX_DCSX_DDNY_H */
