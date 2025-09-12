
#ifndef SRC_WRAPPERS_ROCM_COMMON_ROCSPARSE_H
#define SRC_WRAPPERS_ROCM_COMMON_ROCSPARSE_H

#include <rocsparse/rocsparse.h>

#include "esinfo/eslog.hpp"



#ifndef CHECK
#define CHECK(status) do { _check((status), __FILE__, __LINE__); } while(false)
#endif
inline void _check(rocsparse_status status, const char *file, int line)
{
    if (status != rocsparse_status_success) {
        espreso::eslog::error("rocSPARSE Error %d. In file '%s' on line %d\n", status, file, line);
    }
}



namespace espreso {
namespace gpu {

    template<typename I>
    static inline rocsparse_indextype get_rocsparse_index_type()
    {
        if constexpr(std::is_same_v<I, int32_t>) return rocsparse_indextype_i32;
        if constexpr(std::is_same_v<I, int64_t>) return rocsparse_indextype_i64;
    }

    template<typename T>
    static inline rocsparse_datatype get_rocsparse_data_type()
    {
        if constexpr(std::is_same_v<T, float>)  return rocsparse_datatype_f32_r;
        if constexpr(std::is_same_v<T, double>) return rocsparse_datatype_f64_r;
        if constexpr(std::is_same_v<T, std::complex<float>>)  return rocsparse_datatype_f32_c;
        if constexpr(std::is_same_v<T, std::complex<double>>) return rocsparse_datatype_f64_c;
    }

    static inline rocsparse_order get_rocsparse_order(char order)
    {
        if(order == 'R') return rocsparse_order_row;
        if(order == 'C') return rocsparse_order_column;
        eslog::error("wrong order\n");
    }

}
}



namespace espreso {
namespace gpu {
namespace spblas {

    struct _handle
    {
        rocsparse_handle h;
        hipStream_t get_stream()
        {
            hipStream_t stream;
            CHECK(rocsparse_get_stream(h, &stream));
            return stream;
        }
    };

}
}
}




#endif /* SRC_WRAPPERS_ROCM_COMMON_ROCSPARSE_H */
