
#ifndef SRC_WRAPPERS_CUDA_COMMON_CUBLAS_H
#define SRC_WRAPPERS_CUDA_COMMON_CUBLAS_H

#include <cublas_v2.h>

#include "esinfo/eslog.hpp"



#ifndef CHECK
#define CHECK(status) do { _check((status), __FILE__, __LINE__); } while(false)
#endif
inline void _check(cublasStatus_t status, const char *file, int line)
{
    if (status != CUBLAS_STATUS_SUCCESS) {
        espreso::eslog::error("CUBLAS Error %d %s: %s. In file '%s' on line %d\n", status, cublasGetStatusName(status), cublasGetStatusString(status), file, line);
    }
}



template<typename T> struct cpp_to_cublas_type { using type = T; };
template<> struct cpp_to_cublas_type<std::complex<float>> { using type = cuComplex; };
template<> struct cpp_to_cublas_type<std::complex<double>> { using type = cuDoubleComplex; };
template<typename T> using cpp_to_cublas_type_t = typename cpp_to_cublas_type<T>::type;



namespace espreso {
namespace gpu {
namespace dnblas {

    struct _handle
    {
        cublasHandle_t h;
    };

}
}
}



#endif /* SRC_WRAPPERS_CUDA_COMMON_CUBLAS_H */
