
#ifndef SRC_WRAPPERS_CUDA_OPERATIONS_W_CUSPARSE_GEMM_DCSX_DDNY_DDNZ_H
#define SRC_WRAPPERS_CUDA_OPERATIONS_W_CUSPARSE_GEMM_DCSX_DDNY_DDNZ_H

#include "gpu/operations/gemm_dcsx_ddny_ddnz.h"



namespace espreso {
namespace gpu {
namespace operations {



#ifdef HAVE_CUDA



template<typename T, typename I>
struct w_cusparse_gemm_dcsx_ddny_ddnz_data;

template<typename T, typename I>
class w_cusparse_gemm_dcsx_ddny_ddnz : public gemm_dcsx_ddny_ddnz<T,I>
{
public:
    w_cusparse_gemm_dcsx_ddny_ddnz();
    virtual ~w_cusparse_gemm_dcsx_ddny_ddnz();
protected:
    using gemm_dcsx_ddny_ddnz<T,I>::q;
    using gemm_dcsx_ddny_ddnz<T,I>::handle_spblas;
    using gemm_dcsx_ddny_ddnz<T,I>::A;
    using gemm_dcsx_ddny_ddnz<T,I>::B;
    using gemm_dcsx_ddny_ddnz<T,I>::C;
    using gemm_dcsx_ddny_ddnz<T,I>::ws_persistent;
    using gemm_dcsx_ddny_ddnz<T,I>::wss_internal;
    using gemm_dcsx_ddny_ddnz<T,I>::wss_persistent;
    using gemm_dcsx_ddny_ddnz<T,I>::wss_tmp_preprocess;
    using gemm_dcsx_ddny_ddnz<T,I>::wss_tmp_perform;
    using gemm_dcsx_ddny_ddnz<T,I>::alpha;
    using gemm_dcsx_ddny_ddnz<T,I>::beta;
protected:
    void internal_setup() override;
    void internal_preprocess(void * ws_tmp) override;
    void internal_perform(void * ws_tmp) override;
private:
    std::unique_ptr<w_cusparse_gemm_dcsx_ddny_ddnz_data<T,I>> data;
};



#else



template<typename T, typename I>
class w_cusparse_gemm_dcsx_ddny_ddnz : public gemm_dcsx_ddny_ddnz<T,I>
{
public:
    w_cusparse_gemm_dcsx_ddny_ddnz() { eslog::error("cuda wrapper is not available\n"); }
    virtual ~w_cusparse_gemm_dcsx_ddny_ddnz() = default;
};



#endif



}
}
}

#endif /* SRC_WRAPPERS_CUDA_OPERATIONS_W_CUSPARSE_GEMM_DCSX_DDNY_DDNZ_H */
