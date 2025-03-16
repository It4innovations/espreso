
#ifndef SRC_WRAPPERS_CUDA_OPERATIONS_W_CUDA_SUBMATRIX_DCSX_DCSX_H
#define SRC_WRAPPERS_CUDA_OPERATIONS_W_CUDA_SUBMATRIX_DCSX_DCSX_H

#include "gpu/operations/submatrix_dcsx_dcsx.h"



namespace espreso {
namespace gpu {
namespace operations {



#ifdef HAVE_CUDA



template<typename T, typename I>
class w_cuda_submatrix_dcsx_dcsx : public submatrix_dcsx_dcsx<T,I>
{
public:
    w_cuda_submatrix_dcsx_dcsx() = default;
    virtual ~w_cuda_submatrix_dcsx_dcsx() = default;
protected:
    virtual void internal_setup() override;
    virtual void internal_preprocess(void * ws_tmp) override;
    virtual void internal_perform(void * ws_tmp) override;
private:
    size_t wss_pers_startptrs = 0;
    size_t wss_pers_endptrs = 0;
    size_t wss_pers_outptrs = 0;
    I * src_start_ptrs = nullptr;
    I * src_end_ptrs = nullptr;
    I * dst_ptrs = nullptr;
};



#else



template<typename T, typename I>
class w_cuda_submatrix_dcsx_dcsx : public submatrix_dcsx_dcsx<T,I>
{
public:
    w_cuda_submatrix_dcsx_dcsx() { eslog::error("cuda wrapper is not available\n"); }
    virtual ~w_cuda_submatrix_dcsx_dcsx() = default;
};



#endif



}
}
}

#endif /* SRC_WRAPPERS_CUDA_OPERATIONS_W_CUDA_SUBMATRIX_DCSX_DCSX_H */
