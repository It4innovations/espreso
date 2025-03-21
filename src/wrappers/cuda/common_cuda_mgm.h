
#ifndef SRC_WRAPPERS_CUDA_COMMON_CUDA_MGM_H
#define SRC_WRAPPERS_CUDA_COMMON_CUDA_MGM_H

#include <cuda_runtime.h>
#include <cuComplex.h>
#include <complex>

#include "esinfo/eslog.hpp"



#ifndef CHECK
#define CHECK(status) do { _check((status), __FILE__, __LINE__); } while(false)
#endif
inline void _check(cudaError_t error_code, const char *file, int line)
{
    if (error_code != cudaSuccess)
    {
        espreso::eslog::error("CUDA Error %d %s: %s. In file '%s' on line %d\n", error_code, cudaGetErrorName(error_code), cudaGetErrorString(error_code), file, line);
    }
}



template<typename T> struct cpp_to_cuda_type { using type = T; };
template<> struct cpp_to_cuda_type<std::complex<float>> { using type = cuComplex; };
template<> struct cpp_to_cuda_type<std::complex<double>> { using type = cuDoubleComplex; };
template<typename T> using cpp_to_cuda_type_t = typename cpp_to_cuda_type<T>::type;



namespace espreso {
namespace gpu {
namespace mgm {

    struct _device
    {
        int gpu_idx;
    };

    struct _queue
    {
        cudaStream_t stream;
    };

}
}
}



#endif /* SRC_WRAPPERS_CUDA_COMMON_CUDA_MGM_H */
