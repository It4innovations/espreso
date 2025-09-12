
#ifndef SRC_WRAPPERS_ROCM_OPERATIONS_W_ROCBLAS_TRSM_DDNX_DDNY_H
#define SRC_WRAPPERS_ROCM_OPERATIONS_W_ROCBLAS_TRSM_DDNX_DDNY_H

#include "gpu/operations/trsm_ddnx_ddny.h"



namespace espreso {
namespace gpu {
namespace operations {



#ifdef HAVE_ROCM



struct w_rocblas_trsm_ddnx_ddny_data;

template<typename T>
class w_rocblas_trsm_ddnx_ddny : public trsm_ddnx_ddny<T>
{
public:
    w_rocblas_trsm_ddnx_ddny();
    virtual ~w_rocblas_trsm_ddnx_ddny();
protected:
    using trsm_ddnx_ddny<T>::q;
    using trsm_ddnx_ddny<T>::handle_dnblas;
    using trsm_ddnx_ddny<T>::A;
    using trsm_ddnx_ddny<T>::X;
    using trsm_ddnx_ddny<T>::wss_tmp_perform;
protected:
    void internal_setup() override;
    void internal_perform(void * ws_tmp) override;
private:
    void do_call(char stage);
private:
    std::unique_ptr<w_rocblas_trsm_ddnx_ddny_data> data;
};



#else



template<typename T>
class w_rocblas_trsm_ddnx_ddny : public trsm_ddnx_ddny<T>
{
public:
    w_rocblas_trsm_ddnx_ddny() { eslog::error("rocm wrapper is not available\n"); }
    virtual ~w_rocblas_trsm_ddnx_ddny() = default;
};



#endif



}
}
}

#endif /* SRC_WRAPPERS_ROCM_OPERATIONS_W_ROCBLAS_TRSM_DDNX_DDNY_H */
