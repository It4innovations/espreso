
#ifndef SRC_WRAPPERS_CUDA_OPERATIONS_W_CUBLAS_HERK_DDNX_DDNY_H
#define SRC_WRAPPERS_CUDA_OPERATIONS_W_CUBLAS_HERK_DDNX_DDNY_H

#include "gpu/operations/herk_ddnx_ddny.h"



namespace espreso {
namespace gpu {
namespace operations {



#ifdef HAVE_CUDA



struct w_cublas_herk_ddnx_ddny_data;

template<typename T>
class w_cublas_herk_ddnx_ddny
{
public:
    w_cublas_herk_ddnx_ddny();
    virtual ~w_cublas_herk_ddnx_ddny();
protected:
    virtual void internal_setup() override;
    virtual void internal_perform(void * /*ws_tmp*/) override;
private:
    std::unique_ptr<w_cublas_herk_ddnx_ddny_data> data;
};



#else



template<typename T>
class w_cublas_herk_ddnx_ddny : public herk_ddnx_ddny<T>
{
public:
    w_cublas_herk_ddnx_ddny() { eslog::error("cuda wrapper is not available\n"); }
    virtual ~w_cublas_herk_ddnx_ddny() = default;
};



#endif



}
}
}



#endif /* SRC_WRAPPERS_CUDA_OPERATIONS_W_CUBLAS_HERK_DDNX_DDNY_H */
