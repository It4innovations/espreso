
#ifndef SRC_WRAPPERS_CUDA_W_CUDA_COMMON_H_
#define SRC_WRAPPERS_CUDA_W_CUDA_COMMON_H_

#include <cuda_runtime.h>

#include "esinfo/eslog.h"



#ifndef CHECK
#define CHECK(status) do { _check((status), __FILE__, __LINE__); } while(false)
#endif

inline static void _check(cudaError_t error_code, const char *file, int line)
{
    if (error_code != cudaSuccess)
    {
        char str[1000];
        snprintf(str, sizeof(str), "CUDA Error %d %s: %s. In file '%s' on line %d\n", error_code, cudaGetErrorName(error_code), cudaGetErrorString(error_code), file, line);
        espreso::eslog::error(str);
    }
}

#endif
