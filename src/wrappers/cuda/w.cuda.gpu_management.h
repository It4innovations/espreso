
#ifndef SRC_WRAPPERS_CUDA_W_CUDA_GPU_MANAGEMENT_H_
#define SRC_WRAPPERS_CUDA_W_CUDA_GPU_MANAGEMENT_H_

#include <cuda_runtime.h>

#include "esinfo/eslog.h"



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

#endif /* SRC_WRAPPERS_CUDA_W_CUDA_GPU_MANAGEMENT_H_ */
