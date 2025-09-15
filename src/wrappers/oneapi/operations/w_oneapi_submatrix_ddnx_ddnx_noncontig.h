
#ifndef SRC_WRAPPERS_ONEAPI_OPERATIONS_W_ONEAPI_SUBMATRIX_DDNX_DDNX_H
#define SRC_WRAPPERS_ONEAPI_OPERATIONS_W_ONEAPI_SUBMATRIX_DDNX_DDNX_H

#include "gpu/operations/submatrix_ddnx_ddnx_noncontig.h"



namespace espreso {
namespace gpu {
namespace operations {



#ifdef HAVE_ONEAPI



template<typename T, typename I>
class w_oneapi_submatrix_ddnx_ddnx_noncontig : public submatrix_ddnx_ddnx_noncontig<T,I>
{
public:
    w_oneapi_submatrix_ddnx_ddnx_noncontig() = default;
    virtual ~w_oneapi_submatrix_ddnx_ddnx_noncontig() = default;
protected:
    using submatrix_ddnx_ddnx_noncontig<T,I>::q;
    using submatrix_ddnx_ddnx_noncontig<T,I>::d_M_src;
    using submatrix_ddnx_ddnx_noncontig<T,I>::d_M_dst;
    using submatrix_ddnx_ddnx_noncontig<T,I>::d_row_map;
    using submatrix_ddnx_ddnx_noncontig<T,I>::d_col_map;
protected:
    void internal_perform() override;
};



#else



template<typename T, typename I>
class w_oneapi_submatrix_ddnx_ddnx_noncontig : public submatrix_ddnx_ddnx_noncontig<T,I>
{
public:
    w_oneapi_submatrix_ddnx_ddnx_noncontig() { eslog::error("oneapi wrapper is not available\n"); }
    virtual ~w_oneapi_submatrix_ddnx_ddnx_noncontig() = default;
};



#endif



}
}
}

#endif /* SRC_WRAPPERS_ONEAPI_OPERATIONS_W_ONEAPI_SUBMATRIX_DDNX_DDNX_H */
