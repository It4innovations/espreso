
#ifndef SRC_WRAPPERS_CUDA_OPERATIONS_W_CUBLAS_TRSM_DDNX_DDNY_H
#define SRC_WRAPPERS_CUDA_OPERATIONS_W_CUBLAS_TRSM_DDNX_DDNY_H

#include "gpu/operations/trsm_ddnx_ddny.h"



namespace espreso {
namespace gpu {
namespace operations {



#ifdef HAVE_CUDA



template<typename T>
class w_cublas_trsm_ddnx_ddny : public trsm_ddnx_ddny<T>
{
public:
    w_cublas_trsm_ddnx_ddny();
    virtual ~w_cublas_trsm_ddnx_ddny();
protected:
    using trsm_ddnx_ddny<T>::q;
    using trsm_ddnx_ddny<T>::handle_dnblas;
    using trsm_ddnx_ddny<T>::A;
    using trsm_ddnx_ddny<T>::X;
    using trsm_ddnx_ddny<T>::wss_tmp_perform;
protected:
    void internal_setup() override;
    void internal_perform(void * ws_tmp) override;
};



#else



template<typename T>
class w_cublas_trsm_ddnx_ddny : public trsm_ddnx_ddny<T>
{
public:
    w_cublas_trsm_ddnx_ddny() { eslog::error("cuda wrapper is not available\n"); }
    virtual ~w_cublas_trsm_ddnx_ddny() = default;
};



#endif



}
}
}

#endif /* SRC_WRAPPERS_CUDA_OPERATIONS_W_CUBLAS_TRSM_DDNX_DDNY_H */
