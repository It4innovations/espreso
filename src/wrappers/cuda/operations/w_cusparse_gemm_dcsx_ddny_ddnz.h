
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
    virtual void internal_set_matrix_A() override;
    virtual void internal_set_matrix_B() override;
    virtual void internal_set_matrix_C() override;
    virtual void internal_setup() override;
    virtual void internal_preprocess(void * ws_tmp) override;
    virtual void internal_perform(void * ws_tmp) override;
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
