
#ifndef SRC_WRAPPERS_ROCM_OPERATIONS_W_ROCBLAS_GEMM_DDNX_DDNY_DDNZ_H
#define SRC_WRAPPERS_ROCM_OPERATIONS_W_ROCBLAS_GEMM_DDNX_DDNY_DDNZ_H

#include "gpu/operations/gemm_ddnx_ddny_ddnz.h"



namespace espreso {
namespace gpu {
namespace operations {



#ifdef HAVE_ROCM



struct w_rocblas_gemm_ddnx_ddny_ddnz_data;

template<typename T>
class w_rocblas_gemm_ddnx_ddny_ddnz : public gemm_ddnx_ddny_ddnz<T>
{
public:
    w_rocblas_gemm_ddnx_ddny_ddnz();
    virtual ~w_rocblas_gemm_ddnx_ddny_ddnz();
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
    void internal_setup() override;
    void internal_perform(void * ws_tmp) override;
private:
    void do_call(char stage);
private:
    std::unique_ptr<w_rocblas_gemm_ddnx_ddny_ddnz_data> data;
};



#else



template<typename T>
class w_rocblas_gemm_ddnx_ddny_ddnz : public gemm_ddnx_ddny_ddnz<T>
{
public:
    w_rocblas_gemm_ddnx_ddny_ddnz() { eslog::error("rocm wrapper is not available\n"); }
    virtual ~w_rocblas_gemm_ddnx_ddny_ddnz() = default;
};



#endif



}
}
}

#endif /* SRC_WRAPPERS_ROCM_OPERATIONS_W_ROCBLAS_GEMM_DDNX_DDNY_DDNZ_H */
