
#ifndef SRC_WRAPPERS_ROCM_COMMON_ROCBLAS_H
#define SRC_WRAPPERS_ROCM_COMMON_ROCBLAS_H

#include <rocblas/rocblas.h>

#include "esinfo/eslog.hpp"



#ifndef CHECK
#define CHECK(status) do { _check((status), __FILE__, __LINE__); } while(false)
#endif
inline void _check(rocblas_status status, const char *file, int line)
{
    if (status != rocblas_status_success && status != rocblas_status_size_unchanged && status != rocblas_status_size_increased) {
        espreso::eslog::error("rocBLAS Error %d. In file '%s' on line %d\n", status, file, line);
    }
}



template<typename T> struct cpp_to_rocblas_type { using type = T; };
template<> struct cpp_to_rocblas_type<std::complex<float>> { using type = rocblas_float_complex; };
template<> struct cpp_to_rocblas_type<std::complex<double>> { using type = rocblas_double_complex; };
template<typename T> using cpp_to_rocblas_type_t = typename cpp_to_rocblas_type<T>::type;



namespace espreso {
namespace gpu {
namespace dnblas {

    struct _handle
    {
        rocblas_handle h;
    };

}
}
}



#endif /* SRC_WRAPPERS_ROCM_COMMON_ROCBLAS_H */
