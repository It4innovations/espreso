
#ifndef SRC_WRAPPERS_CUDA_COMMON_CUSPARSE_H
#define SRC_WRAPPERS_CUDA_COMMON_CUSPARSE_H

#include <cusparse.h>

#include "esinfo/eslog.hpp"



#ifndef CHECK
#define CHECK(status) do { _check((status), __FILE__, __LINE__); } while(false)
#endif
inline void _check(cusparseStatus_t status, const char *file, int line)
{
    if (status != CUSPARSE_STATUS_SUCCESS) {
        espreso::eslog::error("CUSPARSE Error %d %s: %s. In file '%s' on line %d\n", status, cusparseGetErrorName(status), cusparseGetErrorString(status), file, line);
    }
}



template<typename I>
inline cusparseIndexType_t cusparse_index_type()
{
    if constexpr(std::is_same_v<I, int32_t>) return CUSPARSE_INDEX_32I;
    if constexpr(std::is_same_v<I, int64_t>) return CUSPARSE_INDEX_64I;
}

template<typename T>
inline cudaDataType_t cusparse_data_type()
{
    if constexpr(std::is_same_v<T, float>)  return CUDA_R_32F;
    if constexpr(std::is_same_v<T, double>) return CUDA_R_64F;
    if constexpr(std::is_same_v<T, std::complex<float>>)  return CUDA_C_32F;
    if constexpr(std::is_same_v<T, std::complex<double>>) return CUDA_C_64F;
}

inline cusparseOrder_t cusparse_order(char order)
{
    if(order == 'R') return CUSPARSE_ORDER_ROW;
    if(order == 'C') return CUSPARSE_ORDER_COL;
    espreso::eslog::error("wrong order\n");
}



namespace espreso {
namespace gpu {
namespace spblas {

    struct _handle
    {
        cusparseHandle_t h;
        cudaStream_t get_stream()
        {
            cudaStream_t stream;
            CHECK(cusparseGetStream(h, &stream));
            return stream;
        }
    };

}
}
}

#endif /* SRC_WRAPPERS_CUDA_COMMON_CUSPARSE_H */
