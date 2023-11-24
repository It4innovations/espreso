
#pragma once

#ifdef MY_ROC



#include <cstdio>

#include <hip/hip_runtime.h>
#include <rocsparse/rocsparse.h>
#include <rocblas/rocblas.h>
#include <rocsolver/rocsolver.h>

#include "my_common.hpp"
#include "matrices.hpp"
#include "cbmba.hpp"

namespace {


#ifndef CHECK
#define CHECK(status) do { check((status), __FILE__, __LINE__); } while(false)
inline static void check(hipError_t error_code, const char *file, int line)
{
    if (error_code != hipSuccess)
    {
        fprintf(stderr, "HIP Error %d %s: %s. In file '%s' on line %d\n", error_code, hipGetErrorName(error_code), hipGetErrorString(error_code), file, line);
        fflush(stderr);
        exit(error_code);
    }
}
inline static void check(rocsparse_status status, const char *file, int line)
{
    if (status != rocsparse_status_success)
    {
        fprintf(stderr, "ROCSPARSE Error %d: %s. In file '%s' on line %d\n", status, "no tostring I could find", file, line);
        fflush(stderr);
        exit(status);
    }
}
inline static void check(rocblas_status status, const char *file, int line)
{
    if (status != rocblas_status_success && status != rocblas_status_size_unchanged && status != rocblas_status_size_increased)
    {
        fprintf(stderr, "ROCBLAS Error %d: %s. In file '%s' on line %d\n", status, "no tostring I could find", file, line);
        fflush(stderr);
        exit(status);
    }
}
// rocsolver_status is deprecated, use rocblas_status
#endif

template<typename I>
rocsparse_indextype rocsparse_index_type()
{
    if constexpr(std::is_same_v<I, int32_t>) return rocsparse_indextype_i32;
    if constexpr(std::is_same_v<I, int64_t>) return rocsparse_indextype_i64;
}

template<typename T>
rocsparse_datatype rocsparse_data_type()
{
    if constexpr(std::is_same_v<T, float>)  return rocsparse_datatype_f32_r;
    if constexpr(std::is_same_v<T, double>) return rocsparse_datatype_f64_r;
}

}

#include "hip_stuff_common.hpp"

#endif
