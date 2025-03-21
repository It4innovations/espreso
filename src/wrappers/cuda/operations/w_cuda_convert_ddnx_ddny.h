
#ifndef SRC_WRAPPERS_CUDA_OPERATIONSW_CUDA_COPY_DDNX_DDNX_H
#define SRC_WRAPPERS_CUDA_OPERATIONSW_CUDA_COPY_DDNX_DDNX_H

#include "gpu/operations/convert_ddnx_ddny.h"



namespace espreso {
namespace gpu {
namespace operations {



#ifdef HAVE_CUDA



template<typename T>
class w_cuda_convert_ddnx_ddny : public convert_ddnx_ddny<T>
{
public:
    w_cuda_convert_ddnx_ddny() = default;
    virtual ~w_cuda_convert_ddnx_ddny() = default;
protected:
    using convert_ddnx_ddny<T>::q;
    using convert_ddnx_ddny<T>::handle_dnblas;
    using convert_ddnx_ddny<T>::M_src;
    using convert_ddnx_ddny<T>::M_dst;
    using convert_ddnx_ddny<T>::wss_tmp_perform;
protected:
    virtual void internal_setup() override;
    virtual void internal_perform(void * ws_tmp) override;
};



#else



template<typename T>
class w_cuda_convert_ddnx_ddny : public submatrix_dcsx_dcsx<T>
{
public:
    w_cuda_convert_ddnx_ddny() { eslog::error("cuda wrapper is not available\n"); }
    virtual ~w_cuda_convert_ddnx_ddny() = default;
};



#endif



}
}
}

#endif /* SRC_WRAPPERS_CUDA_OPERATIONSW_CUDA_COPY_DDNX_DDNX_H */
