
#ifndef SRC_WRAPPERS_ROCM_W_ROCM_GPU_MANAGEMENT_H_
#define SRC_WRAPPERS_ROCM_W_ROCM_GPU_MANAGEMENT_H_

#include <hip/hip_runtime.h>

#include "esinfo/eslog.h"



#ifndef CHECK
#define CHECK(status) do { _check((status), __FILE__, __LINE__); } while(false)
#endif

inline void _check(hipError_t error_code, const char *file, int line)
{
    if (error_code != hipSuccess)
    {
        char str[1000];
        snprintf(str, sizeof(str), "HIP Error %d %s: %s. In file '%s' on line %d\n", error_code, hipGetErrorName(error_code), hipGetErrorString(error_code), file, line);
        espreso::eslog::error(str);
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
        hipStream_t stream;
    };

}
}
}

#endif /* SRC_WRAPPERS_ROCM_W_ROCM_GPU_MANAGEMENT_H_ */
