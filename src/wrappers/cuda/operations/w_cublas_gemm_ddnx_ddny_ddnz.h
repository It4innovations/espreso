
#ifndef SRC_WRAPPERS_CUDA_OPERATIONS_W_CUBLAS_GEMM_DDNX_DDNY_DDNZ_H
#define SRC_WRAPPERS_CUDA_OPERATIONS_W_CUBLAS_GEMM_DDNX_DDNY_DDNZ_H

#include "gpu/operations/gemm_ddnx_ddny_ddnz.h"



namespace espreso {
namespace gpu {
namespace operations {



#ifdef HAVE_CUDA



struct w_cublas_gemm_ddnx_ddny_ddnz_data;

template<typename T>
class w_cublas_gemm_ddnx_ddny_ddnz : public gemm_ddnx_ddny_ddnz<T>
{
public:
    w_cublas_gemm_ddnx_ddny_ddnz();
    virtual ~w_cublas_gemm_ddnx_ddny_ddnz();
protected:
    using gemm_ddnx_ddny_ddnz<T>::q;
    using gemm_ddnx_ddny_ddnz<T>::handle_dnblas;
    using gemm_ddnx_ddny_ddnz<T>::A;
    using gemm_ddnx_ddny_ddnz<T>::B;
    using gemm_ddnx_ddny_ddnz<T>::C;
    using gemm_ddnx_ddny_ddnz<T>::wss_tmp_perform;
    using gemm_ddnx_ddny_ddnz<T>::alpha;
    using gemm_ddnx_ddny_ddnz<T>::beta;
protected:
    virtual void internal_setup() override;
    virtual void internal_perform(void * ws_tmp) override;
private:
    std::unique_ptr<w_cublas_gemm_ddnx_ddny_ddnz_data> data;
};



#else



template<typename T>
class w_cublas_gemm_ddnx_ddny_ddnz : public gemm_ddnx_ddny_ddnz<T>
{
public:
    w_cublas_gemm_ddnx_ddny_ddnz() { eslog::error("cuda wrapper is not available\n"); }
    virtual ~w_cublas_gemm_ddnx_ddny_ddnz() = default;
};



#endif



}
}
}

#endif /* SRC_WRAPPERS_CUDA_OPERATIONS_W_CUBLAS_GEMM_DDNX_DDNY_DDNZ_H */
