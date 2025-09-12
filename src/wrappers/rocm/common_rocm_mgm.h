
#ifndef SRC_WRAPPERS_ROCM_COMMON_ROCM_MGM_H
#define SRC_WRAPPERS_ROCM_COMMON_ROCM_MGM_H

#include <hip/hip_runtime.h>

#include "esinfo/eslog.hpp"



#ifndef CHECK
#define CHECK(status) do { _check((status), __FILE__, __LINE__); } while(false)
#endif
inline void _check(hipError_t error_code, const char *file, int line)
{
    if (error_code != hipSuccess)
    {
        espreso::eslog::error("HIP Error %d %s: %s. In file '%s' on line %d\n", error_code, hipGetErrorName(error_code), hipGetErrorString(error_code), file, line);
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



#endif /* SRC_WRAPPERS_ROCM_COMMON_ROCM_MGM_H */
